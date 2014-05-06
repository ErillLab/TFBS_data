import MySQLdb
import pandas as pd
import numpy as np
import re

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
    return pd.Series(dict(genome_accession=genome_acc, start_pos=start_pos,
                          end_pos=end_pos, strand=strand))

def read_sites(con):
    """Read sites into a data frame."""
    site_df = pd.read_sql("SELECT * FROM sites WHERE func_site_type_guid=668", con)
    # drop unnecessary columns
    drop_columns = ['pkg_guid', 'fl_real_name', 'func_site_type_guid',
                    'struct_site_type_guid', 'regulator_guid',
                    'fl_dna_rna', 'pos_from', 'pos_from_guid',
                    'pfo_type_id', 'pfo_side_guid', 'pos_to',
                    'pos_to_guid', 'pto_type_id', 'pto_side_guid',
                    'site_len', 'signature', 'last_update']
    site_df = site_df.drop(drop_columns, axis=1)

    # parse description field to get genome accession and location
    site_df = pd.merge(site_df.drop('descript', axis=1),
                       site_df.apply(lambda x: parse_site_desc_field(x['descript']), axis=1),
                       left_index=True,
                       right_index=True)
    # parse name field to get the TF name
    
    
    # exclude rows with no genome_acc, site location or strand info
    site_df = site_df.dropna()
    return site_df

def read_genomes(con):
    """Read genomes"""
    genome_df = pd.read_sql("SELECT * FROM dict_genomes", con)
    genome_df = genome_df.rename(columns={'name': 'species'}) # column renaming
    return genome_df

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
    return groups

def mysqldump_to_csv(con):
    """Read database tables to dataframes"""
    art_df = read_articles(con)
    genome_df = read_genomes(con)
    site_df = read_sites(con)
    # merge site_df and art_df to get PMIDs
    site_df = pd.merge(site_df, art_df, on='art_guid')
    # merge site_df and genome_df to get species names
    site_df = pd.merge(site_df, genome_df, on='genome_guid')

    # sort by genome_accession
    site_df = site_df.sort(['genome_accession', 'name'])
    site_df.index = range(1, len(site_df)+1)
    
    # Write to csv
    cols_to_write = ['genome_accession',
                     'name',
                     'start_pos',
                     'end_pos',
                     'strand',
                     'sequence',
                     'pmid',
                     ]

    site_df.to_csv('regtransbase_mysql.csv', sep=',', cols=cols_to_write)

def parse_fasta_seq((desc, seq)):
    """Parse decription line to get genome accession and site location"""
    match = re.match(r">(\S+) \|(\w+):(\d+)-(\d+) \[.+\]", desc)
    if not match:
        print "parse error on"
        print desc
        print seq
        return None

    groups = match.groups()
    return {'TF': groups[0],
            'genome_accession': groups[1],
            'start_pos': groups[2],
            'end_pos': groups[3],
            'sequence': seq,
           }

def fasta_to_csv(fasta_file="regtransbase.fasta"):
    """Given the fasta file containing binding sites, parse the description
    lines and write the sites to the csv file"""
    with open (fasta_file, 'r') as f:
        lines = map(lambda l: l.strip(), f.readlines())
    # parse fasta file
    fasta_seqs = zip(lines[::2], lines[1::2])
    data = [d for d in map(parse_fasta_seq, fasta_seqs) if d]
    df = pd.DataFrame.from_dict(data)
    cols_to_write = ['genome_accession',
                     'TF',
                     'start_pos',
                     'end_pos',
                     'sequence']
    df.to_csv('regtransbase_fasta.csv', sep=',', cols=cols_to_write)
    
