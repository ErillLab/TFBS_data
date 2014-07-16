from __future__ import print_function
import csv,os,re
import sys
sys.path.append('..')
from common import *

data_dir = 'data'
header_fields = "CDS,CDS gene name,CDS gene module,Type,Mutant,Auto,Role,Target Gene,Target gene name,target gene module,Motif known,Evidence,PubmedIDs,Binding motif,Is CDS sigma factor".split(',')
# this is a very redundant way of doing it...
cds_index = header_fields.index("CDS")
tf_index = header_fields.index("CDS gene name")
mode_index = header_fields.index("Role")
gene_index = header_fields.index("Target gene name")
evidence_index = header_fields.index("Evidence")
seq_index = header_fields.index("Binding motif")
pubmedIDs_index =header_fields.index("PubmedIDs")

def interpret_mode(mode):
    if mode == 'A':
        return "activator"
    elif mode == 'R':
        return "repressor"
    else:
        return "undefined"

def main():
    outfile = open("coryneRegNet.tsv",'w')
    print(header.strip(),file=outfile) # get rid of trailing newline
    for entry in os.listdir(data_dir):
        sub_dir = os.path.join(data_dir,entry)
        try:
            os.listdir(sub_dir)
        except OSError:
            print("got OSError trying to list:",sub_dir,"continuing")
            continue # skip over non-directories
        organism_name = entry
        print(organism_name)
        ### get genebank filename, downloading it if necessary
        genome_accession = dl_gbk(organism_name)
        gbk_filename = os.path.join('data',organism_name,genome_accession + ".gbk")
        _,fname = os.path.split(gbk_filename)
        genome_accession,ext = os.path.splitext(fname)
        print("genbank filename:",gbk_filename)
        ### end of genebank concerns
        ### get fna if necessary
        fna_filename = dl_fna(organism_name)
        versioned_accession = versioned_accession_from_fna(fna_filename)
        ### end of fna concerns
        d = {}
        gene_regulations_file = os.path.join(sub_dir,'gene_regulations.gb')
        print("starting on:",gene_regulations_file)
        with open(gene_regulations_file) as f:
            all_lines = list(csv.reader(f,delimiter='\t'))
            print(len(all_lines),"lines")
            head,lines = all_lines[0],all_lines[1:]
            for line in lines:
                #print("processing:",line)
                for field in "cds,tf,mode,gene,evidence,seq,pubmedIDs".split(','):
                    exec_string = "%s = line[%s_index]" % (field,field)
                    #print(exec_string)
                    exec(exec_string)
                db_name = "coryneRegNet"
                if not evidence == "experimental":
                    continue
                if not cds in d:
                    tf_accession = protein_accessions_from_locus_tag(cds,gbk_filename)
                    d[cds] = tf_accession
                tf_accession = d[cds]
                if re.match('[acgtACGT]',seq):
                    print("found seq:",seq)
                    interpreted_mode = interpret_mode(mode)
                    interpreted_gene = gene if gene else '-'
                    for raw_site in seq.split(';'): # line contains multiple sites...
                        site = raw_site.upper()
                        #start_ref,stop_ref,strand_ref = find_site_ref(site,fna_filename)
                        start,stop,strand = find_site(site,fna_filename)
                        if (start,stop,strand) == (None,None,None):
                            continue # ignore unlocatable sites or multiple occurrences
                        #assert start_ref == start and stop_ref == stop and strand_ref == strand
                        ufr,dfr = find_flanking_regions(start,stop,strand,site,fna_filename) if strand else (None,None)
                        #print(versioned_accession,tf,tf_accession,ufr,site,dfr,start,stop,strand,gene,mode,db_name,pubmedIDs,sep=',')
                        #print(versioned_accession,tf,tf_accession,ufr,site,dfr,start,stop,strand,interpreted_gene,interpreted_mode,db_name,pubmedIDs,sep='\t',file=outfile)
                        outline = template.substitute(genome_accession=versioned_accession,
                                          TF=tf,
                                          TF_accession=tf_accession,
                                          site_start=start,#gbk_start,
                                          site_end=stop,#gbk_stop,
                                          site_strand=strand,
                                          left_flanking=ufr,
                                          site_sequence=site,
                                          right_flanking=dfr, 
                                          regulated_operon=interpreted_gene,
                                          mode=interpreted_mode,
                                          evidence=pubmedIDs,
                                          database=db_name,
                                          alternative_database_id="-")
                        print(outline,end="\n",file=outfile)
    outfile.close()
