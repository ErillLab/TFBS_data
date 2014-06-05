"""
Reformat exports from various TFBS databases.

The standard format is as follows:
- Genome accession, with the version number (e.g. NC_000913.2)
- TF name
- TF protein accession number (e.g. NP_415648)
- Left and right flanking regions (100bp)
- Position [start, stop), 0-indexed
- Strand {+, -}
- Site sequence
- Regulated operon (as reported in the original database export)
- Mode {Activator, Repressor, Dual}
- Evidence (e.g. list of used experimental techniques, pubmed ids, etc.)
- Database
- Alternative database id
"""

import os
import pandas as pd
from Bio import SeqIO
from Bio import Entrez
Entrez.email = 'sefa1@umbc.edu'

GENOMES_DIR = '/home/sefa/Downloads/'

def reverse_complement(seq):
    """Reverse complement of a sequence"""
    base_complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T':'A',
                        'a': 't', 'c': 'g', 'g': 'c', 't':'a'}
    complement = ''.join(base_complements[l] for l in seq)
    return complement[::-1]

def download_genome(accession, data_dir):
    """Download genome and save it."""
    net_handle = Entrez.efetch(db='nuccore', id=accession,
                               retmode='gbwithparts', rettype='text')
    out_handle = open(os.path.join(data_dir, accession+'.gb'), 'w')
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()

def read_genome(accession, genome_dir=GENOMES_DIR):
    """Read genome from the file"""
    print "Reading", accession
    genome_file = os.path.join(genome_dir, accession + '.gb')
    if not os.path.isfile(genome_file):
        print 'Genome file is not locally available, downloading from NCBI.'
        download_genome(accession, genome_dir)
    record = SeqIO.read(genome_file, 'gb')
    return record

def get_sequence(genome_seq, start, end, strand):
    """Return the sequence from the genome at [start:end). If it is on the
    reverse strand, return the reverse complement"""
    site_seq = genome_seq[start:end]
    if strand == -1:
        site_seq = reverse_complement(site_seq)
    return site_seq

def from_collectf_loc(start, end, strand):
    """Convert collectf location format to my format (0 indexed)"""
    return (start, end+1)

def from_collectf(row, genomes):
    """Convert every row of data from CollecTF csv to the format described at
    the top of this file."""
    genome = genomes[row.genome_accession]
    start, end, strand = row.site_start, row.site_end, row.site_strand
    my_start, my_end = from_collectf_loc(start, end, strand)
    genome_seq = genome.seq.tostring()
    site_seq = get_sequence(genome_seq, my_start, my_end, strand)
    assert site_seq == row.sequence
    left_flanking = get_sequence(genome_seq, my_start-100, my_start, strand)
    right_flanking = get_sequence(genome_seq, my_end, my_end+100, strand)
    return pd.Series(dict(genome_accession=row.genome_accession,
                          TF=row.TF,
                          TF_accession=row.TF_accession,
                          site_start=my_start,
                          site_end=my_end,
                          site_strand=strand,
                          left_flanking=left_flanking,
                          site_sequence=site_seq,
                          right_flanking=right_flanking,
                          regulated_operon=row['regulated genes (locus_tags)'],
                          mode='',
                          evidence=row.experimental_evidence,
                          database='CollecTF',
                          alternative_database_id=''))

def collectf_reformat(collectf_file):
    """Reformat collectf export"""
    df = pd.read_csv(collectf_file, sep='\t')
    genome_accessions = df.genome_accession.unique()
    genomes = dict((acc, read_genome(acc)) for acc in genome_accessions)
    my_df = df.apply(from_collectf, axis=1, args=(genomes,))
    return my_df

def from_mtbreglist_loc(start, end, strand):
    """Convert MtbRegList location format to my format (0 indexed)"""
    return (start, end)

def from_mtbreglist_mode(mode):
    """Map MtbRegList TF mode to my format"""
    modes = {'Pos': 'activator',
             'Neg': 'repressor',
             'Pos/Neg': 'dual',
             'Undef': 'undefined'}
    return modes[mode]

def from_mtbreglist(row, genome):
    """Convert every row of data from MtbRegList to the format described at the
    top of this file"""
    locus_tag = row['TF locus_tag']
    acc, gene = get_tf_accession_and_gene(genome, locus_tag)
    start, end = row['start'], row['end']
    strand = 1 if row['strand'] == '+' else -1
    genome_seq = genome.seq.tostring()
    my_start, my_end = from_mtbreglist_loc(start, end, strand)
    site_seq = get_sequence(genome_seq, my_start, my_end, strand)
    # There are two rows which don't have sites, sequence for these two are
    # inferred using the location.
    if not pd.isnull(row['site']):
        assert site_seq == row['site']
    left_flanking = get_sequence(genome_seq, my_start-100, my_start, strand)
    right_flanking = get_sequence(genome_seq, my_end, my_end+100, strand)
    return pd.Series(dict(genome_accession=row['Genome'],
                          TF=gene,
                          TF_accession=acc,
                          site_start=my_start,
                          site_end=my_end,
                          site_strand=strand,
                          left_flanking=left_flanking,
                          site_sequence=site_seq,
                          right_flanking=right_flanking,
                          regulated_operon='',
                          mode=from_mtbreglist_mode(row.Trx_Impact),
                          evidence=row.Reference,
                          database='MtbRegList',
                          alternative_database_id=''))

def get_tf_accession_and_gene(genome, locus_tag):
    """Given the genome and a locus tag, return the protein name and its RefSeq
    accession number"""
    cdss = [f for f in genome.features if f.type == 'CDS']
    for cds in cdss:
        if cds.qualifiers['locus_tag'][0] == locus_tag:
            acc = cds.qualifiers['protein_id'][0]
            if 'gene' in cds.qualifiers:
                gene = cds.qualifiers['gene'][0]
            else:
                gene = cds.qualifiers['locus_tag'][0]
            return (acc, gene)

def mtbreglist_reformat(mtbreglist_file):
    """Reformat MtbRegList export"""
    df = pd.read_csv(mtbreglist_file, sep=',')
    genome_accessions = df['Genome'].unique()
    # There is only one genome in this database: NC_000962.3
    assert genome_accessions == 'NC_000962.3'
    genome = read_genome('NC_000962.3')
    my_df = df.apply(from_mtbreglist, axis=1, args=(genome,))
    return my_df

def get_tf_accession(genome, tf_name):
    """Given TF name, return its protein accession number"""
    cdss = [f for f in genome.features if f.type == 'CDS']
    for cds in cdss:
        if ('gene' in cds.qualifiers and
            cds.qualifiers['gene'][0].lower() == tf_name.lower()):
            acc = cds.qualifiers['protein_id'][0]
            return acc
    return '-'

def from_regtransbase(row, genomes):
    """Convert every row of data from RegTransBase to the format described at
    the top of this file"""
    genome = genomes[row.genome_accession]
    tf_accession = get_tf_accession(genome, row.TF)
    print row.genome_accession, row.TF, tf_accession
    return row

def regtransbase_reformat(regtransbase_file):
    """Reformat RegTransBase export"""
    df = pd.read_csv(regtransbase_file, sep=',')
    genome_accessions = df['genome_accession'].unique()
    genomes = dict((acc, read_genome(acc)) for acc in genome_accessions)
    # ignore rows with non-RefSeq genome accessions
    df = df[df.genome_accession.str.startswith('NC_')]
    my_df = df.apply(from_regtransbase, axis=1, args=(genomes,))
    return my_df

def merge_all():
    """Merge all data into one csv."""
    mtbreglist_file = 'mtbreglist/mtbreglist.csv'
    mtbreglist_df = mtbreglist_reformat(mtbreglist_file)

    collectf_file = 'collectf/collectfdb.tsv'
    collectf_df = collectf_reformat(collectf_file)

    df = pd.concat([collectf_df, mtbreglist_df])
    cols= ['genome_accession',
           'TF',
           'site_start',
           'site_end',
           'site_strand',
           'left_flanking',
           'site_sequence',
           'right_flanking',
           'regulated_operon',
           'mode',
           'evidence',
           'database',
           'alternative_database_id']
    df.to_csv('tfbs_data_merged.tsv', cols=cols, sep='\t', index=False)

