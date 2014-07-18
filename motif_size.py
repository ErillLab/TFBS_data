"""Motif size for each TF/genome"""

import pandas as pd

def group_sites(df):
    """Group binding sites by TF accession and genome accession"""
    # There are duplicate entries for some binding sites (i.e. two binding sites
    # on forward and reverse strand, remove duplicates
    df.drop_duplicates(cols=['TF_accession',
                             'genome_accession',
                             'site_start',
                             'site_end'],
                       inplace=True)

    grps = df.groupby(['TF_accession', 'genome_accession', 'TF'])
    return grps

def main():
    df = pd.read_csv('tfbs_data_merged.tsv', sep='\t')
    groups = group_sites(df)
    motif_sizes = groups.size().reset_index().rename(columns={0:'motif_size'})
    motif_sizes.to_csv('motif_sizes.csv', index=False)
    
