"""Motif size for each TF/genome"""

import pandas as pd
import os
from Bio import Entrez
from Bio import SeqIO
Entrez.email = 'sefa1@umbc.edu'
import matplotlib.pyplot as plt
from ggplot import *

def get_org_name(acc):
    """Get genome record"""
    downloads_dir = '/home/sefa/Downloads'
    genbank_file = os.path.join(downloads_dir, acc+'.gb')
    # if not exist, download from NCBI genome database
    if not os.path.isfile(genbank_file):
        print "Downloading genbank file", acc
        net_handle = Entrez.efetch(db='nuccore', id=acc,
                                   retmode='gbwithparts', rettype='text')
        out_handle = open(genbank_file, 'w')
        out_handle.write(net_handle.read())
        net_handle.close()
        out_handle.close()
    # read record
    record = SeqIO.read(genbank_file, 'genbank')
    return record.description


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

def motif_sizes():
    df = pd.read_csv('tfbs_data_merged.tsv', sep='\t')
    groups = group_sites(df)
    motif_sizes = groups.size().reset_index().rename(columns={0:'motif_size'})
    motif_sizes.to_csv('motif_sizes.csv', index=False)
    motif_sizes = motif_sizes[motif_sizes['motif_size'] <= 100]
    p = ggplot(aes(x='motif_size'), motif_sizes) + geom_histogram(binwidth=1)
    ggsave(p, "motif_sizes.pdf")
    


def num_sites_by_db():
    plt.clf()
    df = pd.read_csv('tfbs_data_merged.tsv', sep='\t')
    x = df.groupby('database').size()
    x.plot(kind='bar', rot=0)
    plt.xlabel("")
    plt.ylabel("Number of sites")
    plt.savefig('num_sites_by_db.pdf')
