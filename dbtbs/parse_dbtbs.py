from __future__ import print_function
import xml.etree.ElementTree as ET
import re
import sys
sys.path.append('..')
from common import *

### Preliminaries
fna_filename = 'genome_data/Bacillus_subtilis_168_uid57675/NC_000964.fna'
gbk_filename = 'genome_data/Bacillus_subtilis_168_uid57675/NC_000964.gbk'
### Forget about sigma factors for now
#relevant_fields = "tfac gene sigma regulation sequence".split()
relevant_fields = "tfac gene regulation sequence".split()
#printable_fields ="tfac_text gene_text sigma_text regulation_text site".split()
printable_fields ="tfac_text gene_text regulation_text site".split()
    
def main():
    tree = ET.parse("dbtbs.xml")
    root = tree.getroot()
    csv = open("dbtbs.csv","w")
    tf_accession_dict = {}
    versioned_accession = versioned_accession_from_fna(fna_filename)
    print(header.strip(),sep=',',file=csv)
    for i,promoter in enumerate(root.iter('promoter')):
        for field in relevant_fields:
            exec_string1 = "%s = promoter.find('%s')" % (field,field)
            exec(exec_string1)
            exec_string2 = "%s_text = eval('%s').text if eval('%s') is not None else None" % (field,field,field)
            exec(exec_string2)
        if any([eval("%s_text" % field) is None for field in relevant_fields]):
            print("FAILURE | missing relevant field:",tfac_text,gene_text,regulation_text,sequence_text)
            #print("FAILURE | missing relevant field:",tfac_text,gene_text,sigma_text,regulation_text,sequence_text)
            continue
        sites = [site.replace("/","") for site in re.findall(r"\{(.*?)\}",sequence_text)]
        if not tfac_text in tf_accession_dict:
            print("could not find:",tfac_text,"; searching")
            tf_accession_dict[tfac_text] = protein_accessions_from_gene_name(decap(tfac_text),gbk_filename)
        else:
            print("found",tfac_text," in dictionary")
        tf_accession = tf_accession_dict[tfac_text]
        if tf_accession is None:
            print("FAILURE | Could not find:",tfac_text)
            continue
        tf = tfac_text
        print("gene:",gene)
        gene = gene_text
        print("gene:",gene)
        if regulation_text == "Positive":
            mode = "A"
        elif regulation_text == "Negative":
            mode = "R"
        else:
            print("WARNING | unrecognized regulation mode:",regulation_text)
            mode = regulation_text
        db_name = "tbdbs"
        alternate_db_id = ""
        evidence = "[" + ";".join([elm.text for elm in (promoter.find('reference').iter())][1:]) + "]"
        for site in sites:
            if len(site) < 7:
                print("FAILURE | site too short at:",len(site))
                continue
            start,stop,strand = find_site(site,fna_filename)
            if (start,stop,strand) == (None,None,None):
                print("FAILURE | multiple or no matches")
                continue
            ufr,dfr = find_flanking_regions(start,stop,strand,site,fna_filename) if strand else (None,None)
            print("SUCCESS | ",versioned_accession,tf,tf_accession,ufr,site,dfr,start,stop,
                  strand,gene,mode,db_name,evidence,sep=',')
            print(versioned_accession,tf,tf_accession,ufr,site,dfr,start,stop,
                  strand,gene,mode,db_name,evidence,alternate_db_id,sep=',',file=csv)
            #print(",".join(eval(field_name) for field_name in printable_fields),file=csv)
    csv.close()
