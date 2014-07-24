from __future__ import print_function
import csv,re
from utils import get_ecoli_genome,verbose_gen,wc
from parse_tf_set import gene_name_from_tf
from Bio import SeqIO
import sys
sys.path.append('..')
from common import *

gbk_fname = "NC_000913.gbk"
fna_fname = "NC_000913.2.fna"
gbk = SeqIO.read(gbk_fname,'genbank')

with open("BindingSiteSet.txt") as f:
    lines = list(csv.reader(f,delimiter='\t'))

with open("NC_000913.2.fna") as f:
    ecoli_v2 = "".join(map(lambda line:line.strip(),f.readlines())[1:])
rev_ecoli_v2 = wc(ecoli_v2)
G = len(ecoli_v2)

db_fields = "tf_id,tf_name,bs_id,left_pos,right_pos,strand,interaction_id,regulandum,mode,promoter,tss_pos,site,evidence".split(',')

for field in db_fields:
    exec_string = "%s=%s" % (field,db_fields.index(field))
    exec(exec_string)

def from_my_index(my_start_pos,site,strand):
    if strand == "forward":
        my_stop_pos = my_start_pos + len(site)
        their_start_pos,their_stop_pos = my_start_pos + 1, my_stop_pos
    elif strand == "reverse":
        my_stop_pos = my_start_pos + len(upper_site)
        their_start_pos,their_stop_pos = my_start_pos-len(upper_site)+1,my_start_pos
    return their_start_pos,their_stop_pos

def from_their_index(their_start,their_stop,strand):
    """Convert a position in regulonDB indexing to mine (more important function).
    Function satisfies following relations:

    if strand == 'forward':
        my_start,my_stop = from_their_index(their_start,their_stop,strand)
        site == ecoli[my_start:my_stop]
    if strand == 'reverse':
        my_start,my_stop = from_their_index(their_start,their_stop,strand)
        site == rev_ecoli[my_start:my_stop] == ecoli[G - my_start:G - my_stop]
    """
    if strand == "forward":
        my_start , my_stop = their_start - 1,their_stop
    elif strand == "reverse":
        #rev_ecoli_v2[G-stop:G-start+1] #OLD
        #rev_ecoli_v2[G-their_stop:G-their_start+1] #Mon Jul 14 21:53:21 EDT 2014
        their_stop = their_start + len(upper_site)
        my_start,my_stop = G - their_stop + 1,G - their_start + 1
    else:
        print("didn't recognize strand:",strand)
    return my_start,my_stop

def to_gbk_index(my_start,my_stop,strand):
    if strand == "forward":
        return my_start + 1,my_stop + 1
    elif strand == "reverse":
        return my_start,my_stop

def interpret_mode(mode):
    mode_dict =  {"+":"activator",
                  "-":"repressor",
                  "+-":"dual"}
    if mode in mode_dict:
        return mode_dict[mode]
    else:
        return "undefined"

def main():
    tf_dict = {}
    with open("regulonDB.tsv",'w') as outfile:
        outfile.write(header + "\n")
        for line in verbose_gen(lines,modulus=100):
            #print line
            if line[0].startswith('#'):
                continue
            raw_site = line[site]
            upper_site = filter(lambda c:c.isupper(),raw_site)
            if not upper_site:
                continue
            matches = find_site(upper_site,fna_fname,return_all=True)
            if len(matches) == 1:
                start,stop,site_strand = matches[0]
            elif len(matches) > 1:
                raw_matches = find_site(raw_site.upper(),fna_fname,return_all=True)
                if len(raw_matches) > 1:
                    print("WARNING | too many raw matches:",raw_matches)
                    print("compare to:",start,stop)
                    raw_matches = [sorted(raw_matches,key=lambda(x,y,z):abs(x-start))[0]]
                    print("chose:",raw_matches)
                    # fall through to next if block since len(raw_matches) == 1 now...
                if len(raw_matches) == 1:
                    raw_start,raw_stop,raw_strand = raw_matches[0]
                    start,stop,site_strand = find(lambda (x,y,z): (raw_start <= x and
                                                                   y <= raw_stop and
                                                                   z == raw_strand),matches)
                    print("found unique match by raw matching")
                else:
                    print("FAILURE | no raw matches")
                    continue
            else:
                print("FAILURE | no matches")
                continue
                
            ufr,dfr = find_flanking_regions(start,stop,site_strand,upper_site,fna_fname) if site_strand else (None,None)
            tf = line[tf_name]
            if not tf in tf_dict:
                tf_dict[tf] = protein_accessions_from_gene_name(tf,gbk_fname)
                print("assigned TF accession: %s -> %s" % (tf,tf_dict[tf]))
            tf_accession = tf_dict[tf]
            outline = template.substitute(genome_accession="NC_000913.2", # Build v. 2
                                          TF = line[tf_name],
                                          TF_accession=tf_accession,
                                          site_start=start,#gbk_start,
                                          site_end=stop,#gbk_stop,
                                          site_strand=site_strand,
                                          left_flanking=ufr,
                                          site_sequence=upper_site,
                                          right_flanking=dfr, 
                                          regulated_operon=line[regulandum],
                                          mode=interpret_mode(line[mode]),
                                          evidence=line[evidence].replace(',',';'),
                                          database="regulonDB",
                                          alternative_database_id=line[bs_id])
            print("SUCCESS")
            print(outline,end="\n",file=outfile)
    print("finished")
