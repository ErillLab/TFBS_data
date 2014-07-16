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
        print "didn't recognize strand:",strand 
    return my_start,my_stop

def to_gbk_index(my_start,my_stop,strand):
    if strand == "forward":
        return my_start + 1,my_stop + 1
    elif strand == "reverse":
        return my_start,my_stop

def main():
    tf_dict = {}
    with open("regulonDB.csv",'w') as outfile:
        outfile.write(header)
        for line in verbose_gen(lines,modulus=100):
            #print line
            if line[0].startswith('#'):
                continue
            raw_site = line[site]
            upper_site = filter(lambda c:c.isupper(),raw_site)
            if not upper_site:
                continue
            # start = int(line[left_pos])
            # stop = int(line[right_pos])
            # my_start,my_stop = from_their_index(start,stop,line[strand])
            # gbk_start,gbk_stop = to_gbk_index(my_start,my_stop,line[strand])
            # if line[strand] == 'forward':
            #     genome = ecoli_v2
            #     site_strand = 1
            # elif line[strand] == 'reverse':
            #     genome = rev_ecoli_v2
            #     site_strand = -1
            # assert upper_site == genome[my_start:my_stop]
            # ufr = genome[my_start-100:my_start]
            # dfr = genome[my_stop:my_stop+100]
            #assert ufr+upper_site+dfr in genome (works,but slow)
            matches = find_site(upper_site,fna_fname,return_all=True)
            if len(matches) == 1:
                start,stop,site_strand = matches[0]
            elif len(matches) > 1:
                raw_matches = find_site(raw_site.upper(),fna_fname,return_all=True)
                if len(raw_matches) > 1:
                    print "WARNING | too many raw matches:",raw_matches
                    print "compare to:",start,stop
                    raw_matches = [sorted(raw_matches,key=lambda(x,y,z):abs(x-start))[0]]
                    print "chose:",raw_matches
                    # fall through to next if block since len(raw_matches) == 1 now...
                if len(raw_matches) == 1:
                    raw_start,raw_stop,raw_strand = raw_matches[0]
                    start,stop,site_strand = find(lambda (x,y,z): (raw_start <= x and
                                                                   y <= raw_stop and
                                                                   z == raw_strand),matches)
                    print "found unique match by raw matching"
                else:
                    print "FAILURE | no raw matches"
                    continue
            else:
                print("FAILURE | no matches")
                continue
                
            ufr,dfr = find_flanking_regions(start,stop,site_strand,upper_site,fna_fname) if site_strand else (None,None)
            if not tf_name in tf_dict:
                tf_dict[tf_name] = protein_accessions_from_gene_name(decap(line[tf_name]),gbk_fname)
            tf_accession = tf_dict[tf_name]
            outline = template.substitute(genome_accession="NC_000913.2", # Build v. 2
                                          tf_name = line[tf_name],
                                          tf_accession=tf_accession,
                                          ufr=ufr,
                                          site=upper_site,
                                          dfr=dfr,
                                          start_pos=start,#gbk_start,
                                          stop_pos=stop,#gbk_stop,
                                          strand=site_strand,
                                          operon=line[regulandum],
                                          mode=line[mode],
                                          db_name="regulonDB",
                                          evidence=line[evidence],
                                          alternate_db_id=line[bs_id])
            print "SUCCESS"
            outfile.write(outline)    
    print "finished"
