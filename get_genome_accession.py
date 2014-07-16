#!/usr/bin/python2.6

# This script was cannabalized from the dl_genomes script originally
# written for the reciprocal_blast project.

import ftputil
import string
import os
import sys
sys.path.append('..')
from utils import *

#Where to put the genomes (script will create sub directories
#using the same names as the NCBI use).  This directory must
#exist already:
base_path="."

host = ftputil.FTPHost('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')
host.chdir('/genomes/Bacteria/')

def org_matches_dir(org, org_dir):
    """Returns whether org_dir is the org_dir of org"""
    return all(word.lower() in org_dir.lower() for word in org.split('_'))
    
def genome_accession_from_species_name(species_name): 
    print "Searching for " + species_name
    ### First check to see if we've downloaded it already
    print "checking data directory for gbk file"
    data_dir_contents = os.listdir('data')
    subdir = find(lambda fname:fname.startswith(species_name),data_dir_contents)
    if not subdir is None:
        print "Found sub-directory:",subdir
        subdir_contents = os.listdir(os.path.join('data',subdir))
        gbk_file = find(lambda fname:fname.endswith("gbk"),subdir_contents)
        if not gbk_file is None:
            print "Found gbk file:",gbk_file
            accession, ext = os.path.splitext(gbk_file)
            return accession
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
    if len(unique_accessions) == 1:
        return accessions[0]
    else:
        print "Found multiple accesions"
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
        return unique_accessions[int(choice)]
    
