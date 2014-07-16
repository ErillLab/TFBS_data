from Bio import SeqIO
from Bio import Entrez
Entrez.email = "pon2@umbc.edu"
from string import Template
from get_genome_accession import genome_accession_from_species_name
from utils import levenshtein,wc,head
import os,re
from utils import find
import ftputil
host = ftputil.FTPHost('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')
host.chdir('/genomes/Bacteria/')

template = Template("$genome_accession, $tf_name, $tf_accession, $ufr, $site, $dfr, $start_pos, $stop_pos, $operon, $mode, $db_name, $evidence, $alternate_db_id\n")
header = "genome_accession, TF, TF_accession, left_flanking, site_sequence, right_flanking, site_start, site_stop, site_strand, regulated_operon, mode, database, evidence, alternative_database_id"


def dl_gbk(species_name):
    """Return genomic accession for species, downloading gbk file if necessary"""
    print "searching for gbk for:",species_name
    data_dir_contents = os.listdir('data')
    subdir = species_name
    super_dir = os.path.join('data',subdir)
    if not subdir is None:
        print "Found sub-directory:",subdir
        subdir_contents = os.listdir(os.path.join('data',subdir))
        gbk_file = find(lambda fname:fname.endswith("gbk"),subdir_contents)
        if not gbk_file is None:
            print "Found gbk file:",gbk_file
            accession,ext = os.path.splitext(gbk_file)
            print "accession:",accession
            return accession
        else:
            print "didn't find gbk file in:",subdir
            print "attempting to download"
    host.chdir('/genomes/Bacteria/')
    dir_list = host.listdir(host.curdir)
    inserted_list = sorted(dir_list + [species_name])
    idx = inserted_list.index(species_name)
    best_guesses = inserted_list[idx-10:idx+10]
    for i,dir_name in enumerate(dir_list):
        if dir_name in best_guesses:
            print i,dir_name
    best_guess = (list(enumerate(dir_list))[idx])
    print best_guess
    print "choose a directory or [Enter] for best guess: %s,%s" % best_guess
    choice = raw_input()
    if choice == '':
        choice = str(idx)
    dir_name = dir_list[int(choice)]
    print "Interpreted %s as %s" % (species_name,dir_name)
    host.chdir('/genomes/Bacteria/' + dir_name + '/')
    file_list = host.listdir(host.curdir)
    accessions = [f[:f.index('.')] for f in file_list]
    unique_accessions = list(set(accessions))
    if len(unique_accessions) > 1:
        print "Found multiple accessions"
        # Compare size of fnas to choose largest
        filesize = lambda f: host.path.getsize(f)
        sized_fnas = sorted([(filesize(f),f) for f in file_list if f.endswith(".fna")],reverse=True)
        print ".fna filesizes:"
        for fs,f in sized_fnas:
            print fs,f
        largest_size,largest_fna = sized_fnas[0]
        largest_accession = largest_fna[:largest_fna.index('.')]
        print "Choose an accession or [Enter] for largest:",largest_accession
        for i,v in enumerate(unique_accessions):
            print i,v
        choice = raw_input()
        if choice == '':
            choice = str(unique_accessions.index(largest_accession))
        accession = unique_accessions[int(choice)]
    else:
        accession = accessions[0]
    gbk_filename = accession + ".gbk"
    target_path = os.path.join('data',subdir,gbk_filename)
    print "attempting to download:", gbk_filename, "to:",target_path
    host.download(gbk_filename,target_path)
    return accession

def dl_fna(species_name):
    """Dl fna if necessary, return filename"""
    accession = dl_gbk(species_name)
    print "accession:",accession
    fna_name = accession + ".fna"
    print "fna_name:",fna_name
    target_path = os.path.join("data",species_name,fna_name)
    if os.path.isfile(target_path):
        print "found fna:",target_path
        return target_path
    print "didn't find fna for:",species_name,"downloading"
    host.chdir('/genomes/Bacteria/')
    dir_list = host.listdir(host.curdir)
    sorted_dir_list = sorted(dir_list,key=lambda fname:levenshtein(species_name,fname))
    for dir_name in sorted_dir_list:
        print "trying:",dir_name
        try:
            host.chdir('/genomes/Bacteria/' + dir_name + '/')
            sub_dir_list = host.listdir(host.curdir)
            if find(lambda name:name.startswith(accession),sub_dir_list):
                host.download(fna_name,target_path)
                return target_path
        except:
            "something got fucked up!"
            continue
    print "Couldn't find fna for:",species_name
    return None
    
def protein_accessions_from_gene_name(gene_name,gbk_fname):
    """Given a gene name and a gbk filename, return the protein accession"""
    print "Searching for:",gene_name,"in:",gbk_fname
    gbk = SeqIO.read(gbk_fname,'genbank')
    p_accs = []
    for feature in gbk.features:
        if 'gene' in feature.qualifiers and 'protein_id' in feature.qualifiers:
            if gene_name in feature.qualifiers['gene']:
                p_accs.append(feature.qualifiers['protein_id'][0])
    if not p_accs:
        print "couldn't find:",gene_name
        return None
    if len(p_accs) > 1:
        print "Warning, found %s entries for %s: %s" % (len(p_accs),gene_name,p_accs)
    return p_accs[0]

def protein_accessions_from_locus_tag(locus_tag,gbk_fname):
    """Given a locus_tag and a gbk filename, return the protein accession"""
    print "Searching for:",locus_tag,"in:",gbk_fname
    gbk = SeqIO.read(gbk_fname,'genbank')
    p_accs = []
    for feature in gbk.features:
        if 'locus_tag' in feature.qualifiers and 'protein_id' in feature.qualifiers:
            if locus_tag in feature.qualifiers['locus_tag']:
                p_accs.append(feature.qualifiers['protein_id'][0])
    if not p_accs:
        print "couldn't find:",locus_tag
        return None
    if len(p_accs) > 1:
        print "Warning, found %s entries for %s: %s" % (len(p_accs),gene_name,p_accs)
    return p_accs[0]

def genome_from_fna(fna_filename):
    with open(fna_filename) as f:
        genome = "".join([line.strip() for line in f.readlines()[1:]])
    return genome
        
def find_site_ref(site,fna_filename):
    """WRONG: rev matches are indexed backwards"""
    genome = genome_from_fna(fna_filename)
    regexp = re.compile(site)
    fwd_matches = [(lambda (start,stop):(start,stop,+1))(m.span()) for m in regexp.finditer(genome)]
    rev_matches = [(lambda (start,stop):(start,stop,-1))(m.span()) for m in regexp.finditer(wc(genome))]
    matches = fwd_matches + rev_matches
    print matches
    if len(matches) == 1:
        print "found unique match for %s in %s" % (site,fna_filename)
        return head(matches)
    elif len(matches) > 1:
        print "found multiple matches for %s in %s" % (site,fna_filename)
        return head(matches)
    else:
        print "couldn't find' match for %s in %s" % (site,fna_filename)
        return (None,None,None)

def find_site(site,fna_filename,return_all=False):
    genome = genome_from_fna(fna_filename)
    fwd_regexp = re.compile(site)
    rev_regexp = re.compile(wc(site))
    fwd_matches = [(lambda (start,stop):(start,stop,+1))(m.span()) for m in fwd_regexp.finditer(genome)]
    rev_matches = [(lambda (start,stop):(start,stop,-1))(m.span()) for m in rev_regexp.finditer(genome)]
    matches = fwd_matches + rev_matches
    #print matches
    if len(matches) == 1:
        print "found unique match for %s in %s" % (site,fna_filename)
        return head(matches) if not return_all else matches
    elif len(matches) > 1:
        print "found multiple matches for %s in %s" % (site,fna_filename)
        return (None,None,None) if not return_all else matches
    else:
        print "couldn't find' match for %s in %s" % (site,fna_filename)
        return (None,None,None) if not return_all else []

def find_flanking_regions(start,stop,strand,site,fna_filename):
    genome = genome_from_fna(fna_filename)
    offset = 100
    ufr = genome[start-offset:start]
    dfr = genome[stop:stop+offset]
    if strand == -1:
        ufr,dfr = wc(dfr),wc(ufr)
        assert wc(ufr + site + dfr) in genome
    else:
        assert ufr + site + dfr in genome
    return ufr,dfr

def versioned_accession_from_fna(fna_filename):
    with open(fna_filename) as f:
        header = f.readline()
    fields = header.split("|")
    return fields[3]

def decap(xs):
    return xs[0].lower() + xs[1:]
