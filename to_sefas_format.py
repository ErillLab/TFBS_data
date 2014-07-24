"""
Convert Pat's csv's to tsvs.

Sayeth Sefa:
-----
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
-----
"""

import pandas as pd

sefa_cols= ['genome_accession',
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

regulonDB_file = "regulonDB/regulonDB.csv"
coryneRegNet_file = "coryneRegNet/coryneRegNet.csv"
dbtbs_file = "dbtbs/dbtbs.csv"


