"""
This module is used to generate csv files from RegTransBase mysqldump and TFBS
alignments, available at RegTransBase website.

To generate csv file from mysqldump,

  con = connect_db(host, db, user, passwd)
  mysqldump_to_csv(con)

To generate csv file from TFBS alignments

  fasta_to_csv(alignment_dir)

"""

import MySQLdb
import pandas as pd
import re
import os

def connect_db(host="localhost", db="dbregulation_update", user="sefa", passwd=""):
    """Connect mysql database and return a connection"""
    con = MySQLdb.connect(host=host, db=db, user=user, passwd=passwd)
    return con

def read_articles(con):
    """Read articles into a data frame"""
    df = pd.read_sql("SELECT * from articles", con)
    art_df = df.loc[:, ['art_guid', 'pmid']]
    return art_df

def parse_site_desc_field(desc):
    """Parse site description field to extract genome accession number, site
    location and strand."""
    match = re.match(r"NCBI-ACCESSION:(\w+):(\d+):(\d+):(\S+)", desc)
    genome_acc = None
    start_pos = None
    end_pos = None
    strand = None
    if match is None:
        pass
    else:
        genome_acc, start_pos, end_pos, strand = match.groups()
        start_pos = int(start_pos)
        end_pos = int(end_pos)
        start_pos, end_pos = min(start_pos, end_pos), max(start_pos, end_pos)
    return pd.Series(dict(genome_accession=genome_acc, start_pos=start_pos,
                          end_pos=end_pos, strand=strand))

def read_sites(con):
    """Read sites into a data frame."""
    site_df = pd.read_sql("SELECT * FROM sites WHERE func_site_type_guid=668", con)
    # drop unnecessary columns
    drop_columns = ['pkg_guid', 'fl_real_name', 'func_site_type_guid',
                    'struct_site_type_guid',
                    'fl_dna_rna', 'pos_from', 'pos_from_guid',
                    'pfo_type_id', 'pfo_side_guid', 'pos_to',
                    'pos_to_guid', 'pto_type_id', 'pto_side_guid',
                    'site_len', 'signature', 'last_update']
    site_df = site_df.drop(drop_columns, axis=1)
    site_df.rename(columns={'name': 'description'}, inplace=True)     # rename columns
    site_df.sequence = site_df.sequence.apply(lambda x: x.upper())
    # parse description field to get genome accession and location
    site_df = pd.merge(site_df.drop('descript', axis=1),
                       site_df.apply(lambda x: parse_site_desc_field(x['descript']), axis=1),
                       left_index=True,
                       right_index=True)
    # exclude rows with no genome_acc, site location or strand info
    site_df = site_df.dropna()
    # site start and end positions are integers
    site_df.start_pos = site_df.start_pos.astype(int)
    site_df.end_pos = site_df.end_pos.astype(int)
    return site_df

def read_genomes(con):
    """Read genomes"""
    genome_df = pd.read_sql("SELECT * FROM dict_genomes", con)
    genome_df = genome_df.rename(columns={'name': 'species'}) # column renaming
    return genome_df

def read_regulators(con):
    """Read protein names"""
    regulator_df = pd.read_sql("SELECT regulator_guid, name FROM regulators", con)
    regulator_df = regulator_df.rename(columns={'name': 'TF'})
    return regulator_df

"""
In Regtransbase database, sites are linked to articles, not
experiments. However, an article may be linked to many sites, as well as many
experiments. There is no way to be sure that which experiments are for which
sites. Therefore, the experiment information is not merged into the sites table.
"""

def read_experiments(con):
    """Read experiments into a data frame"""
    # read experiments
    exp_df = pd.read_sql("SELECT * FROM experiments", con)
    # The database also contains types for experiments where each experiment may
    # have more than one type.
    # read experiment types
    exp_type_df = pd.read_sql("SELECT * FROM exp2technique_types", con)
    # read experiment names
    exp_name_df = pd.read_sql("SELECT * FROM dict_exp_technique_types", con)
    # join exp types and their names
    exp_type_df = pd.merge(exp_type_df, exp_name_df, on='exp_technique_type_guid')
    # join experiment table and types
    exp_df = pd.merge(exp_df, exp_type_df, on='exp_guid')
    # drop unnecessary columns
    exp_df = exp_df.drop(['pkg_guid', 'last_update'], axis=1)
    return exp_df

def mysqldump_to_csv(con):
    """Read database tables to dataframes"""
    art_df = read_articles(con)
    genome_df = read_genomes(con)
    regulator_df = read_regulators(con)
    site_df = read_sites(con)
    # merge site_df and art_df to get PMIDs
    site_df = pd.merge(site_df, art_df, on='art_guid')
    # merge site_df and genome_df to get species names
    site_df = pd.merge(site_df, genome_df, on='genome_guid')
    # merge site_df and regulator_df to get TF names
    site_df = pd.merge(site_df, regulator_df, on='regulator_guid')

    # sort by genome_accession
    site_df = site_df.sort(['genome_accession', 'TF'])
    site_df.index = range(1, len(site_df)+1)

    # Write to csv
    cols_to_write = ['genome_accession',
                     'TF',
                     'start_pos',
                     'end_pos',
                     'strand',
                     'sequence',
                     'description',
                     'pmid']

    site_df.to_csv('regtransbase_mysql.csv', sep=',', cols=cols_to_write, index=False)
    return site_df

def parse_fasta_seq((desc, seq)):
    """Parse decription line to get genome accession and site location"""
    match = re.match(r">\S+ \|(\w+):(\d+)-(\d+) \[.+\]", desc)
    if not match:
        print "parse error on"
        print desc
        print seq
        return None

    groups = match.groups()
    genome_acc = groups[0]
    start_pos = int(groups[1])
    end_pos = int(groups[2])
    strand = 1
    if start_pos > end_pos:
        start_pos, end_pos = end_pos, start_pos
        strand = -1

    return {'genome_accession': genome_acc,
            'start_pos': start_pos,
            'end_pos': end_pos,
            'strand': strand,
            'sequence': seq.upper(),
           }

def fasta_to_csv(fasta_dir):
    """Given the directory containing fasta files of binding site alignments,
    parse the description lines and write the sites to the csv file"""
    fasta_files = [f for f in os.listdir(fasta_dir) if f.endswith('.fasta')]
    dfs = [] # data frames
    for fasta_file in fasta_files:
        with open(os.path.join(fasta_dir, fasta_file), 'r') as f:
            lines = map(lambda l: l.strip(), f.readlines())
        # parse fasta file
        fasta_seqs = zip(lines[::2], lines[1::2])
        data = [d for d in map(parse_fasta_seq, fasta_seqs) if d]
        df = pd.DataFrame.from_dict(data)
        df['TF'] = fasta_file.split('_')[0]
        dfs.append(df)

    concatenated = pd.concat(dfs)

    #start and end positions are integers
    concatenated.start_pos = concatenated.start_pos.astype(int)
    concatenated.end_pos = concatenated.end_pos.astype(int)
    concatenated.strand = concatenated.strand.astype(int)

    cols_to_write = ['genome_accession',
                     'TF',
                     'start_pos',
                     'end_pos',
                     'strand',
                     'sequence']
    concatenated.to_csv('regtransbase_fasta.csv', sep=',', cols=cols_to_write, index=False)
