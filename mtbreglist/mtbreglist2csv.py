import os
import re
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
Entrez.email = "sefa1@umbc.edu"

MTB_GENOME_ACCESSION = "NC_000962.3"

def head(xs):
    return xs[0]

def get_genome(acc=MTB_GENOME_ACCESSION):
    """Get genome record"""
    genbank_file = acc + '.gb'
    # if not exist, download from NCBI genome database
    if not os.path.isfile(genbank_file):
        net_handle = Entrez.efetch(db='nuccore', id=acc,
                                   retmode='gbwithparts', rettype='text')
        out_handle = open(genbank_file, 'w')
        out_handle.write(net_handle.read())
        net_handle.close()
        out_handle.close()
    # read record
    record = SeqIO.read(genbank_file, 'genbank')
    return record

def get_gene(genome_record, locus_tag):
    """Find the sequence feature for a particular gene"""
    genes = [f for f in genome_record.features
             if f.type == 'gene' and f.qualifiers['locus_tag'][0] == locus_tag]
    assert len(genes) == 1, "Multiple seqfeature for locus_tag %s" % locus_tag
    return head(genes)

def get_prot(genome_record, locus_tag):
    """Find protein feature in Genbank file for a given locus tag"""
    prots = [f for f in genome_record.features
             if f.type == 'CDS' and f.qualifiers['locus_tag'][0] == locus_tag]
    assert len(prots) == 1
    return head(prots)

def repair_mtb_locus_tags(mtb_export_file):
    """Mycobacterium tuberculosis has locus tags in format (Rv\d+). However, for
    some reason, some of the locus tags in MtbRegList export don't match with
    any gene in the NCBI genome record. It seems that some locus tags (RNA
    coding genes) have been updated at some point in NCBI RefSeq database. These
    locus tags in MtbRegList have been replaced with the one in the genbank file
    if start and end positions match.

    The genomic locations for the genes with MTB-starting locus tags can be
    found at: http://tuberculist.epfl.ch/quicksearch.php?gene+name=MTB000016

    """
    replacements = [("MTB000022", "Rvnt19"),
                    ("MTB000034", "Rvnt30"),
                    ("MTB000004", "Rvnt04"),
                    ("MTB000019", "Rvnr01"),
                    ("MTB000032", "Rvnt28"),
                    ("MTB000016", "Rvnt16")]

    with open(mtb_export_file) as f:
        all_text = f.read()
    # Replace old locus tags with new ones
    for (old, new) in replacements:
        all_text = all_text.replace(old, new)
    with open(mtb_export_file, 'w') as f:
        f.write(all_text)

def reverse_complement(seq):
    """Reverse complement of a sequence"""
    return Seq(seq).reverse_complement().tostring()

def sequence_match(exact_seq, ambiguous_seq):
    """Check if the given exact sequence and ambiguous one agree"""
    return all(eb == ab if ab in "ACTG" else True
               for (eb, ab) in zip(exact_seq, ambiguous_seq))

def binding_site_loc(genome_record, fga_locus_tag, fgb_locus_tag,
                     dist5, dist3):
    """Given flanking gene(s), 5' and 3' distances to these genes, find the
    location of the binding site """
    # Get the location of the genes
    ga = get_gene(genome_record, fga_locus_tag)
    loca = (ga.location.start.position, ga.location.end.position)
    gb = get_gene(genome_record, fgb_locus_tag)
    locb = (gb.location.start.position, gb.location.end.position)
    binding_site_start = loca[1] + dist5
    binding_site_end = locb[0] - dist3
    return (binding_site_start, binding_site_end)

def binding_site_loc_one(genome_record, fg_locus_tag, dist5, dist3):
    """Most binding sites are in intergenic region but some can be inside the
    gene. For such cases, MtbRegList database reports only one gene. Based on
    the location of the gene and distances to its 5' and 3' end, find the
    location of the binding site"""
    # Get the location of the gene
    g = get_gene(genome_record, fg_locus_tag)
    loc = (g.location.start.position, g.location.end.position)
    binding_site_start = loc[0] + dist5
    binding_site_end = loc[1] - dist3
    return (binding_site_start, binding_site_end)

def parse_signature(signature):
    """Parse the signature as defined in MtbRegList export and return actual
    sequence where spacer is denoted with sequence of N's.

    The signature has the pattern S@m[/SP/S@m] where S stands for sequence box,
    SP for spacing and m for mismatch count relative to its root pattern. For
    more info, see MtbRegList paper.
    """
    pattern = r'([ACTG]+)@\d+(/(\d+)/([ACTG]+)@\d+)*'
    m = re.match(pattern, signature)
    seqa, _, spacer_len, seqb = m.groups()

    seq = seqa
    if spacer_len and seqb:
        seq += ('N'*int(spacer_len) + seqb)
    return seq

def identify_binding_site(genome_record, df_row):
    """Given a row of the dataframe, by using Genomic Rrgion, Dist 5' and Dist
    3' columns, extract the binding site from the genome sequence"""
    dist_5 = int(df_row["Dist 5'"])
    dist_3 = int(df_row["Dist 3'"])
    strand = df_row["Strand"]
    assert strand in ['+', '-'], "Invalid strand"
    ambiguous_seq = parse_signature(df_row['Signature'])
    m = re.match(r"ir[PQDT]_(\w+)_(\w+)", df_row['Genomic region'])
    if m: # contains two locus tags and the binding site is between these two
          # genes
        locus_tag_a, locus_tag_b = m.groups()
        start, end = binding_site_loc(genome_record, locus_tag_a, locus_tag_b,
                                      dist_5, dist_3)
    else: # contains only one locus tag
        locus_tag = df_row['Genomic region']
        start, end = binding_site_loc_one(genome_record, locus_tag,
                                          dist_5, dist_3)
    # Extract the binding site
    site = genome_record.seq[start:end].tostring()
    if strand == '-':
        site = reverse_complement(site)
    # Make sure that extracted site is correct
    assert sequence_match(site, ambiguous_seq), "Sequences don't match"
    print start, end, strand, site
    return pd.Series({'start': start, 'end': end, 'strand': strand, 'site':site})

def get_protein_accession(genome_record, df_row):
    """Given a row of the data frame, find the protein accession for the TF
    using its locus tag"""
    locus_tag = df_row["Recognized by"]
    prot = get_prot(genome_record, locus_tag)
    prot_acc = prot.qualifiers['protein_id'][0].split('.')[0]
    return pd.Series({'TF accession': prot_acc})

def read_mtbreglist(export_file="results.txt"):
    """Given the tab delimited text file exported from MtbRegList database, read
    the data that contains promoters, transcription factor binding sites and
    terminators."""
    # fix locus tags
    repair_mtb_locus_tags(export_file)
    # read MtbRegList results into a dataframe
    df = pd.read_csv(export_file, sep='\t', skiprows=2)
    # Remove records other than TF binding sites
    df = df[(df['Recognized by'] != "None") & (df['Recognized by'] != "Unknown")]
    # Parse MtbRegList signature and create a new column of binding sites which
    # may contain N's.
    df['ambiguous_seq'] = df.apply(lambda x:parse_signature(x['Signature']),axis=1)
    # Get genome
    genome_rec = get_genome(MTB_GENOME_ACCESSION)

    df.apply(lambda x: get_protein_accession(genome_rec, x), axis=1)
    
    # identify binding sites
    df = pd.merge(df,
                  df.apply(lambda x: identify_binding_site(genome_rec, x), axis=1),
                  left_index=True,
                  right_index=True)

    df = pd.merge(df,
                  df.apply(lambda x: get_protein_accession(genome_rec, x), axis=1),
                  left_index=True,
                  right_index=True)
    
    # rename some columns
    df.rename(columns={'Recognized by': 'TF locus_tag'}, inplace=True)
    df['Genome'] = MTB_GENOME_ACCESSION
    cols_to_write = ['Genome',
                     'TF accession',
                     'TF locus_tag',
                     'Synonyms',
                     'start',
                     'end',
                     'strand',
                     'site',
                     'Trx_Impact',
                     'Reference']
    df.to_csv("mtbreglist.csv", sep=",", cols=cols_to_write, index=False)
    return df

