"""All following imports necessary for script operations"""

from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
import re
import copy
import time
import pdb
import os
import sys
import pickle
os.chdir("/home/nnard1616/Documents/Python26/TOMM Finder/") #local file in which scripts and relevant files are located


"""Following Class necessary to parse Protein Data"""
class Protein(object):
    def __init__(self, bcd0, gid, length, bits, orgQuery, bcd1="", aa_seq="", cterm=0, cXXc=0, gbfile="", org = "", contig = "", start = "", end = "", DNA = "", rev = ""): #bcd --> is it a sagb, c or d protein?
        self.bcd0 = bcd0 # query protein's functional identity
        self.bcd1 = bcd1 # actual protein's functional identity
        self.gid = gid
        self.length = length
        self.bits = bits # bit score
        self.orgQuery = orgQuery #organism of query sequence
        self.aa_seq = aa_seq 
        self.cterm = cterm #count of prolines at c terminal end.  This is relevant for D protein
        self.cXXc = cXXc #count of cXXc motifs.  This is relevant for C protein
        self.gbfile = gbfile #genbank file of protein
        self.org = org #protein's organism
        self.contig = contig #DNA location
        self.start = start #startingDNA location
        self.end = end #ending DNA location
        self.DNA = DNA #dna sequence
        self.rev = rev #denotes whether gene is read in reverse of database sequence
        
"""PREPARATION OF BLAST RESULT DICTIONARIES"""

f = open("orgDictD.txt", "r")
orgDictD = pickle.load(f)
f.close()
f = open("orgDictC.txt", "r")
orgDictC = pickle.load(f)
f.close()
f = open("orgDictB.txt", "r")
orgDictB = pickle.load(f)
f.close()

"""Helper Functions for my Main function"""

def get_gb_file(prot):
        Entrez.email = "nnard1616@gmail.com"
        try:
            handle = Entrez.efetch(db="protein", id=prot.gid, rettype = "gb")
            gbfile = handle.readlines()
        except:
            print "something went wrong for " + prot.org + "'s " + str(prot.gid)
            return prot
        prot.gbfile = gbfile
        return

"""MAIN FUNCTION"""
n = 0
for d in [orgDictD, orgDictC, orgDictB]:
    n += 1
    for org in d:
        print org, n
        for prot in d[org]:
            get_gb_file(prot)
            
