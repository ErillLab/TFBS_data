"""
RegTransBase has two modules:
- a database of regulatory interactions based on literature
- expertly curated database of transcription factor binding sites

These two modules are provided as mysqldump and fasta files, respectively.

The aim is to distinguish predicted sites and experimentally validated binding
fragments.
"""

import pandas as pd

fasta_csv = "regtransbase_fasta.csv"
mysql_csv = "regtransbase_mysql.csv"

def experimental_or_predicted(text):
    """Given the description field of a site instance, determine if it is
    experimentally verified or predicted"""

    text = text.lower()
    experimental_keywords = ["retarded",
                             "protected",
                             "footprint",
                             "bound by",
                             "shifted",
                             "binding fragment",
                             # same words with typos :/
                             "retaded", "retarted", "retrarded",
                             "reatrded", "retareded", "retatrded",
                             "footrprint", "footoprint",
                             ]

    predicted_keywords = ["predicted",
                          "binding site",
                          "bindng site",
                          "binding motif",
                          "regulatory region",
                          "regulatory site",
                          "box",
                          "segment",
                          "binding region",
                          "gene upstream fragment",
                        ]

    if any(kw in text for kw in experimental_keywords):
        return 1
    elif any(kw in text for kw in predicted_keywords):
        return 2
    else:
        #print "not descriptive description: " + text
        return 0

def filter_exp_verified_sites(df):
    """Given the collection of binding sites/regions from mysqldump, identify
    which ones have experimental evidence"""
    f = df.apply(lambda x: experimental_or_predicted(x.description) == 1, axis=1)
    return df[f]

def filter_predicted_sites(df):
    """Given the collection of binding sites/regions from mysqldump, identify
    which ones have experimental evidence"""
    f = df.apply(lambda x: experimental_or_predicted(x.description) == 2, axis=1)
    return df[f]

def find_experimental_evidence(row, exp_verified_df):
    """Check if the site (given as a row from the alignments csv file) has
    experimental evidence in the other data frame."""
    edf = exp_verified_df # shorter name
    edf = edf[(edf['genome_accession'] == row['genome_accession']) &
              (edf['start_pos'] <= row['start_pos']) &
              (edf['end_pos'] >= row['end_pos']) &
              (edf['TF'] == row['TF'])]

    if edf.empty:
        row['pmid'] = None
    else:
        evidence = edf.iloc[0]
        row['pmid'] = evidence['pmid']
    return row

def get_exp_verified_sites_only():
    """
    RegTransBase provides two different types of data: experimentally validated
    binding sites (or regions) as a mysqldump file and manually curated binding
    sites as a fasta file. Since we are interested in experimentally verified
    binding sites only, this function returns a subset of manually curateed
    binding sites if the site or larger region containing the site has
    experimental evidence in data set."""
    # fasta file contains manually curated binding sites
    fasta_df = pd.read_csv(fasta_csv)
    assert all(fasta_df.start_pos <= fasta_df.end_pos), "fasta_df"

    # mysqldump file contains both experimentally verified regions
    # and predicted sites
    mysql_df = pd.read_csv(mysql_csv)
    assert all(mysql_df.start_pos <= mysql_df.end_pos), "mysql_df"
    
    # we want experimentally verified regions only
    exp_mysql_df = filter_exp_verified_sites(mysql_df)
    pred_mysql_df = filter_predicted_sites(mysql_df)

    # Check all binding sites in the manually curated set and export binding
    # sites that have experimental evidence
    #
    sites_with_evidence = \
        fasta_df.apply(lambda x: find_experimental_evidence(x, exp_mysql_df), axis=1).dropna()
    # Check all predicted sites in mysqldump data set and export binding sites
    # that have experimental evidence
    sites_with_evidence_2 = \
        pred_mysql_df.apply(lambda x: find_experimental_evidence(x, exp_mysql_df), axis=1).dropna()
    cols = ['genome_accession', 'TF', 'start_pos', 'end_pos', 'strand',
            'sequence', 'description', 'pmid']
    concat_df = pd.concat([sites_with_evidence, sites_with_evidence_2])
    concat_df = concat_df[cols]
    concat_df.sort(cols, inplace=True)
    concat_df.to_csv("sites_with_exp_evidence.csv", index=False)
    return concat_df

