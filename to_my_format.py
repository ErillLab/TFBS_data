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
from pyutils import list_utils # github.com/sefakilic/pyutils

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


def from_collectf_mode(mode):
    """Map CollecTF TF mode to my format"""
    if pd.isnull(mode):
        return 'undefined'
        
    modes = {'ACT': 'activator',
             'REP': 'repressor',
             'DUAL': 'dual',
             'N/A': 'undefined'}
    return modes[mode]

def from_collectf_loc(start, end, strand):
    """Convert collectf location format to my format (0 indexed)"""
    return (start-1, end)

def from_collectf(row, genomes):
    """Convert every row of data from CollecTF csv to the format described at
    the top of this file."""
    genome = genomes[row.genome_accession]
    start, end, strand = row.site_start, row.site_end, row.site_strand
    my_start, my_end = from_collectf_loc(start, end, strand)
    print start, end, my_start, my_end
    genome_seq = genome.seq.tostring()
    site_seq = get_sequence(genome_seq, my_start, my_end, strand)
    print row.sequence, site_seq, strand
    assert site_seq == row.sequence, 'site location error'
    left_flanking = get_sequence(genome_seq, my_start-100, my_start, strand)
    right_flanking = get_sequence(genome_seq, my_end, my_end+100, strand)
    if strand == -1:
        left_flanking, right_flanking = right_flanking, left_flanking
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
                          mode=from_collectf_mode(row['mode']),
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
    if strand == -1:
        left_flanking, right_flanking = right_flanking, left_flanking
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

    regulondb_file = 'regulonDB/regulonDB.tsv'
    regulondb_df = pd.read_csv(regulondb_file, sep='\t')

    dbtbs_file = 'dbtbs/dbtbs.tsv'
    dbtbs_df = pd.read_csv(dbtbs_file, sep='\t')

    df = pd.concat([collectf_df,
                    mtbreglist_df,
                    regulondb_df,
                    dbtbs_df])

    # Some TF accesssions have version number (e.g. NP_389668.1) and some don't
    # (e.g. NP_389668). Make them have the same format.
    df['TF_accession'] = df.apply(lambda x: x['TF_accession'].split('.')[0],
                                  axis=1)

    # remove duplicates
    df = remove_duplicates(df)

    # sort rows by TF, TF accesion, genome accession and start position
    df = df.sort(['TF', 'TF_accession', 'genome_accession', 'site_start'])

    cols= ['genome_accession',
           'TF',
           'TF_accession',
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
    return df

def remove_duplicates(df):
    """The same site may occur in different databases, sometimes in slightly
    different location (i.e. (x, x+10) vs. (x-1, x+1), etc.). This function
    removes duplicates from the data frame.

    To merge sites, the method used in CollecTF is adopted. If two sites (from
    same TF and same genome, they are merged if the overlap between two sites
    is >75% of the combined site length.
    """
    def overlap_test(rowa, rowb):
        """Given two rows (i.e. binding sites) from the data frame, check if
        they belong to the same TF/species and if so, check if they overlap.
        """
        def get_overlap(loca, locb):
            """Given two locations, return the overlap ratio."""
            overlap_len = max(0, min(loca[1], locb[1]) - max(loca[0], locb[0]))
            return float(overlap_len) / (loca[1]-loca[0]+1)
        loca = (rowa['site_start'], rowa['site_end'])
        locb = (rowb['site_start'], rowb['site_end'])
        ret = (rowa['genome_accession'] == rowb['genome_accession'] and
               rowa['TF'] == rowb['TF'] and
               (get_overlap(loca, locb) + get_overlap(locb, loca))/2 >= 0.75)

        if ret:
            print rowa['genome_accession'], rowa['TF'], loca, locb,
            print rowa['database'], rowb['database']
        return ret

    rows = df.T.to_dict().values()
    unique_rows = list_utils.nub_by(overlap_test, rows)
    return pd.DataFrame(unique_rows)

if __name__ == '__main__':
    merge_all()
