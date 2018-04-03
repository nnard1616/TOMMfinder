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
import pdb
date = ''
for item in time.gmtime()[0:3]:
  if len(str(item)) ==1:
    date += '0'
    date += str(item)
  else:
    date += str(item)
    
directory = "./%s/" % date  #directory in which dictionaries will be saved (named after the date).  Dictionaries may not be necessary if we start using the computational cluster.


try:
  os.chdir("20110804") #moves pointer to that new directory.
except:
  os.mkdir(directory) #creates that directory if it doesn't already exist.
  os.chdir(directory)# moves there after creating it
missingBCD = ['Bacillus thuringiensis serovar israelensis ATCC 35646', 'Microcystis aeruginosa NIES-843', 'Pseudomonas putida KT2440', 'Streptomyces hygroscopicus ATCC 53653', 'Streptomyces hygroscopicus subsp. jinggangensis'] #Reasons that these are missing are that either the NCBI's record is still under construction or the proteins' genes are not proximal to each other.
missingCD = ['Burkholderia pseudomallei 406e', 'Haloarcula marismortui ATCC 43049', 'Halorubrum lacusprofundi ATCC 49239', 'Natronomonas pharaonis DSM 2160']

#those missing in my BCD list

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
        self.rev = rev #denotes whether gene is complement of gene of database sequence
        
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

BCDorgs = open("../BCDorgs.txt", "r").readlines() #use [:-2] to get rid of "\r\n" at the end of each org name
CDorgs = open("../CDorgs.txt", "r").readlines() #use [:-2] to get rid of "\r\n" at the end of each org name
ALLorgs = BCDorgs + CDorgs

"""Helper Functions for my Main function"""

def File_writer(filename, protList):
    f = open(filename + ".txt", "w")
    for prot in protList:
        f.write(">" + "gi" + str(prot.gid) + "|" + prot.org + "|" + prot.bcd0 + "|" + prot.bcd1 + "|" + "\n")
        f.write(prot.aa_seq + "\n")
        f.write("\n")
    f.close()

def org_file_writer(org, clusters):
    if org.find("\\") != -1 or org.find("/") != -1:
        org = org.replace("\\", "_")
        org = org.replace("/", "_")
    f = open(org + ".txt", "w")
    for prots in clusters:
        num = clusters.index(prots) + 1
        for prot in prots:
            f.write(">" + "gi|" + str(prot.gid) + "|" + prot.org + "|" + "Cluster " + str(num) + "|" + prot.bcd0 + "|" + prot.bcd1 + "|" + "\n")
            f.write(prot.aa_seq + "\n")
            f.write("\n")
    f.close()
    return

def get_gb_file(prot):
        Entrez.email = "nnard1616@gmail.com"
        try:
          handle = Entrez.efetch(db="protein", id=prot.gid, rettype = "gb")
        except:
          print "something is wrong, gbfile of " + str(prot.gid)
          get_gb_file(prot)
        gbfile = handle.readlines()
        prot.gbfile = gbfile
        return prot
      
def get_DNA_sequence(prot):
    Entrez.email = "nnard1616@gmail.com"
    if prot.contig == "":
      return
    start = prot.start
    end = prot.end
    try:
      handle = Entrez.efetch(db="nucleotide", id=prot.contig, rettype = "fasta")
    except:
      print "Something is wrong, DNA of " + prot.gid
      get_DNA_sequence(prot)
    seq_record = SeqIO.read(handle, "fasta")
    contigDNA = str(seq_record.seq)
    if prot.start.find("<") != -1 or prot.start.find(">") != -1:
        start = prot.start[1:]
    if prot.end.find("<") != -1 or prot.end.find(">") != -1:
        end = prot.end[1:]
    DNA = contigDNA[int(start)-1:int(end)]
    prot.DNA = DNA
    return prot

def remove_dups(List): #checks for duplicates in a list, and returns a list that has duplicates removed
    List2 = []
    for item in List:
        try:
            List2.index(item)
            print "wow"
        except:
            List2.append(item)
    return List2

def find_location(prot): #finds and assigns the contig 
    gbfile = prot.gbfile
    line = False
    for lin in gbfile:#SEARCH THROUGH THE GB FILE TO FIND THE CHROMOSOMAL/GENOMIC INFORMATION LOCATION WHERE THIS CERTAIN PROTEIN IS FOUND
        if lin.find("coded_by") != -1:
            line = lin
    if line:
        if line.find('coded_by="complement(') != -1: #for complement
            acc, loc = line.split(":")#ACC IS THE CHROMOSOMAL/GENOMIC ACCESSION NUMBER.  EACH PROTEIN MUST HAVE THE SAME CHROMOSOMAL/GENOMIC ACCESSION NUMBER IN ORDER TO TOMM CLUSTERS
            acc = acc[acc.find('(')+1:]
            start, end = loc[0:-3].split("..")
            prot.contig = acc
            prot.start = start
            prot.end = end
            prot.rev = "Complement"
            return
        if line.find('coded_by="join(') != -1: #generally means it is contained in a eukaryota
            return prot.gid, prot.org
        if line.find('coded_by="') != -1: #for non-complement
            acc, loc = line.split(":")
            acc = acc[acc.find('"')+1:]
            start, end = loc[0:-2].split("..")
            prot.contig = acc
            prot.start = start
            prot.end = end
            return
    return prot.gid, prot.org, "uh-oh"

def check_fusions(prot):
        numP = prot.aa_seq[-60:].count("P") #use for counting how many P are in last 60 AA
        cxxc = re.findall("C[A-Z][A-Z]C", prot.aa_seq) #use for counting number of cxxc motifs
        prot.cterm = numP
        prot.cXXc = len(cxxc)
        gbfile = prot.gbfile
        start = False
        end = False
        for line in gbfile:
                if line.find("FEATURES") != -1:
                        start = gbfile.index(line)
                if line.find("ORIGIN") != -1:
                        end = gbfile.index(line)
        features = ""
        if start == False or end == False:
            return
        for line in gbfile[start:end]:
                features += line
        if (features.find("YcaO") !=-1 or features.find("SagD") != -1 or features.find("mcbD") != -1) and not(features.find("cyclo_dehy") != -1 or features.find("mcbB") != -1 or features.find("SagC") != -1) and not(features.find("mcbC") != -1 or features.find("SagB") != -1 or features.find("FMN binding") != -1 or features.find("Nitro_FMN_reductase") != -1 or features.find("putative FMN binding site") != -1): #checks features if it is a D protein
                prot.bcd1 = "D"
        if not(features.find("YcaO") !=-1 or features.find("SagD") != -1 or features.find("mcbD") != -1) and (features.find("cyclo_dehy") != -1 or features.find("mcbB") != -1 or features.find("SagC") != -1) and not(features.find("mcbC") != -1  or features.find("SagB") != -1 or features.find("FMN binding") != -1 or features.find("Nitro_FMN_reductase") != -1 or features.find("putative FMN binding site") != -1):# checks for C protein
                prot.bcd1 = "C"
        if not(features.find("YcaO") !=-1 or features.find("SagD") != -1 or features.find("mcbD") != -1) and not(features.find("cyclo_dehy") != -1 or features.find("mcbB") != -1 or features.find("SagC") != -1) and (features.find("mcbC") != -1  or features.find("SagB") != -1 or features.find("FMN binding") != -1 or features.find("Nitro_FMN_reductase") != -1 or features.find("putative FMN binding site") != -1):# checks for B protein
                prot.bcd1 = "B"
        if (features.find("YcaO") !=-1 or features.find("SagD") != -1 or features.find("mcbD") != -1) and (features.find("cyclo_dehy") != -1 or features.find("mcbB") != -1 or features.find("SagC") != -1) and not(features.find("mcbC") != -1  or features.find("SagB") != -1 or features.find("FMN binding") != -1 or features.find("Nitro_FMN_reductase") != -1 or features.find("putative FMN binding site") != -1):# checks for CD fusion protein
                prot.bcd1 = "CD"
        if not(features.find("YcaO") !=-1 or features.find("SagD") != -1 or features.find("mcbD") != -1) and (features.find("cyclo_dehy") != -1 or features.find("mcbB") != -1 or features.find("SagC") != -1) and (features.find("mcbC") != -1  or features.find("SagB") != -1 or features.find("FMN binding") != -1 or features.find("Nitro_FMN_reductase") != -1 or features.find("putative FMN binding site") != -1):# checks for BC fusion protein
                prot.bcd1 = "BC"
        if (features.find("YcaO") !=-1 or features.find("SagD") != -1 or features.find("mcbD") != -1) and not(features.find("cyclo_dehy") != -1 or features.find("mcbB") != -1 or features.find("SagC") != -1) and (features.find("mcbC") != -1  or features.find("SagB") != -1 or features.find("FMN binding") != -1 or features.find("Nitro_FMN_reductase") != -1 or features.find("putative FMN binding site") != -1):# checks for BD fusion protein
                prot.bcd1 = "BD"
        if (features.find("YcaO") !=-1 or features.find("SagD") != -1 or features.find("mcbD") != -1) and (features.find("cyclo_dehy") != -1 or features.find("mcbB") != -1 or features.find("SagC") != -1) and (features.find("mcbC") != -1  or features.find("SagB") != -1 or features.find("FMN binding") != -1 or features.find("Nitro_FMN_reductase") != -1 or features.find("putative FMN binding site") != -1):# checks for BCD fusion protein
                prot.bcd1 = "BCD"
        




        if prot.bcd1 == "D" and len(cxxc) >= 1 and prot.length >= 450:
            prot.bcd1 = "CD"
        if prot.bcd1 == "C" and (prot.aa_seq[-10:].count("P") >=3):
            prot.bcd1 = "CD"

        if prot.bcd1 == "B" and len(cxxc) > 1 and prot.length >= 450:
            prot.bcd1 = "?BC"
        
        if not(features.find("YcaO") !=-1 or features.find("SagD") != -1 or features.find("mcbD") != -1) and not(features.find("cyclo_dehy") != -1 or features.find("mcbB") != -1 or features.find("SagC") != -1) and not(features.find("mcbC") != -1 or features.find("FMN binding") != -1) and prot.bcd1 == "":# doesn't have any identifier in the gb file
            prot.bcd1 += "?"
            

def group_by_contig( lisOfProts, orgContigDict):
    contigDict = dict()
    for lis in lisOfProts:
        for prot in lis:
            org = prot.org
##            pdb.set_trace()
            try:
                if contigDict[prot.contig]: #is there a key for this contig? yes
                    try:
                        contigDict[prot.contig].index([prot]) #is this protein already in this contig's value list?
                    except:
                        contigDict[prot.contig].append([prot])#add protein to this contig's value list, as it is not in there
            except: #there does not exist a key for this contig
                contigDict[prot.contig] = [[prot]]
    orgContigDict[org] = contigDict
    return orgContigDict
                
    


def clustering (n, lisProtsOfContig):
    if n == 0:
        return lisProtsOfContig
    else:
        n = 0
        for grp1 in tuple(lisProtsOfContig):
            for grp2 in tuple(lisProtsOfContig):
                if grp1 != grp2:
                    comb = grp1 + grp2
                    lisLoc = []

                    for prot in comb:

                        try:#checks for <1 occurrences
                            lisLoc.append(int(prot.start))
                        except:
                            try:
                                lisLoc.append(int(prot.start[1:]))
                            except:
                                continue
                    try:#sometimes there is no location given, as with plasmids
                        minSt = min(lisLoc)
                        maxSt = max(lisLoc)
                    except:
                        return#not worth working with, so return
                    if abs(minSt - maxSt) < 10000: #distance threshold
                        try:
                            lisProtsOfContig.remove(grp1)
                            lisProtsOfContig.remove(grp2)
                        except:
                            continue
                        lisProtsOfContig.append(comb)
                        n = 1
        clustering(n, lisProtsOfContig)

def consolidate_TOMM(prots):
    for prot1 in prots:
        for prot2 in prots:
            if prot1 != prot2 and prot1.gid == prot2.gid and prot1.bcd1 == "CD" and prot2.bcd1 == "CD":
                prots.remove(prot1)
            if prot1 != prot2 and prot1.gid == prot2.gid and prot1.bcd1 == "D" and prot2.bcd1 == "D" and (prot1.bcd0 == "D" or prot2.bcd0 == "D"):
                prots.remove(prot1)
            if prot1 != prot2 and prot1.gid == prot2.gid and prot1.bcd1 == "C" and prot2.bcd1 == "C" and (prot1.bcd0 == "C" or prot2.bcd0 == "C"):
                prots.remove(prot1)

def eliminate_TOMMBCD(clust, TOMMorgsCD):
    
    for prots in tuple(clust):
        grp = []
        for prot in prots:
            if prot.bcd1.find("B") != -1:
                grp.append("B")
            if prot.bcd1.find("C") != -1:
                grp.append("C")
            if prot.bcd1.find("D") != -1:
                grp.append("D")
            if prot.bcd0 == "B" and prot.bcd1 == "?":
                grp.append("B") #this could screw things up
            if prot.bcd0 == "C" and prot.bcd1 == "?":
                grp.append("C")
            if prot.bcd0 == "D" and prot.bcd1 == "?":
                grp.append("D")
        grp1 = []
        for item in grp:
            try:
                grp1.index(item)
            except:
                grp1.append(item)
        grp1.sort()

        if grp1 != ["B", "C", "D"]:
            if grp1 == ["C","D"]:
                try:
                    TOMMorgsCD[prot.org]
                    TOMMorgsCD[prot.org].append(prots)
                except:
                    TOMMorgsCD[prot.org] = [prots]
            clust.remove(prots)

def eliminate_TOMMCD(clust):
     
    for prots in tuple(clust):
        grp = []
        for prot in prots:
            if prot.bcd1.find("C") != -1:
                grp.append("C")
            if prot.bcd1.find("D") != -1:
                grp.append("D")
            if prot.bcd0 == "C" and prot.bcd1 == "?":
                grp.append("C")
            if prot.bcd0 == "D" and prot.bcd1 == "?":
                grp.append("D")
        grp1 = []
        for item in grp:
            try:
                grp1.index(item)
            except:
                grp1.append(item)
        grp1.sort()

        if grp1 != ["C", "D"]:
            clust.remove(prots)

def check_B(prots):
    bProts = []
    
    for prot in prots:
        if prot.bcd0 == "B":
            bProts.append(prot)
    tbProts = tuple(bProts)
##    if prot.org =="uncultured Prochloron sp. 06037A":
##        pdb.set_trace()
    if len(bProts) > 1:
        for prot in tuple(bProts):
            if prot.bcd1 == "?":
                prots.remove(prot)
                bProts.remove(prot)
            elif prot.bcd1 != "?BC" and prot.bcd1 != "B":
                prots.remove(prot)
                bProts.remove(prot)
                
    else:
        for prot in bProts:
            if prot.bcd1 != "?BC" and prot.bcd1 != "B" and prot.bcd1 != "?":
                prots.remove(prot)
                bProts.remove(prot)
    
    if len(bProts) == 0: #means that all B proteins had "?", therefore we must pick the one with highest bit score
        king = 0
        for prot in tbProts:
            if int(prot.bits) > king:
                king = int(prot.bits)
        for prot in tbProts:
            if int(prot.bits) == king:
                prots.append(prot)
    
print "cool"

for d in [orgDictD, orgDictC, orgDictB]:
  for org in d:
##    print org
    if org == "":
        continue
    for prot in d[org]:
      #get_gb_file(prot)#comment out if gbfile and DNA information is already saved.
      find_location(prot)
      check_fusions(prot)
      #get_DNA_sequence(prot) #comment out if gbfile and DNA information is already saved.
##  f = open("temp/orgDict" + prot.bcd0 + ".txt", "w") #for saving the orgDicts after finding genbank and dna information.  Without this, it would be a pain to find teh data everytime.  Once its saved, you can comment this out as well as teh lines that retrieves teh gb and dna info.
##  pickle.dump(d, f)
##  f.close()


"""comment out next lines to 'beans' if gbfile and dna information are already saved in the dictionaries"""


print "beans"



BCD = dict()
CD = dict() #start here


###groups proteins by contig of organism###
for org in orgDictD: #for bcd tomms
    try:
        dees = orgDictD[org] #for an organism
        cees = orgDictC[org]
        bees = orgDictB[org]
        group_by_contig([dees, cees, bees], BCD)
    except:
        continue


    
for org in orgDictD: #for cd tomms
    try:
        BCD[org]
    except:
        try:
            dees = orgDictD[org]
            cees = orgDictC[org]
            group_by_contig([dees, cees], CD)
        except:
            continue

    
###groups proteins in a contig by proximity###
for d in [BCD, CD]:
    for org in d:
        for con in d[org]:
            clustering(1, d[org][con])
        

###used for counting tomms and finding which organisms contain them.###
TOMMcount = 0
TOMMorgsBCD = dict()
TOMMorgsCD = dict()
for org in BCD:
    for con in BCD[org]:
        for prots in BCD[org][con]:
            clust = []

    
            for prot in prots:
                if prot.bcd0 == "B" :
                    clust.append("B")
                if prot.bcd0 == "C":
                    clust.append("C")
                if prot.bcd0 == "D":
                    clust.append("D")
                if prot.bcd1 == "B" :
                    clust.append("B")
                if prot.bcd1 == "C":
                    clust.append("C")
                if prot.bcd1 == "D":
                    clust.append("D")
                if prot.bcd1 == "CD":
                    clust.append("C")
                    clust.append("D")
            clustNoDups = []

            for let in clust:
                try:
                    clustNoDups.index(let)
                except:
                    clustNoDups.append(let)
            clustNoDups.sort()
            
            if clustNoDups == ["B", "C", "D"]:
                try:
                    TOMMorgsBCD[org]
                    
                    TOMMorgsBCD[org].append(prots)
                except:
                    TOMMorgsBCD[org] = [prots]
            if clustNoDups == ["C", "D"]:
                try:
                    TOMMorgsCD[org]
                    TOMMorgsCD[org].append(prots)
                except:
                    TOMMorgsCD[org] = [prots]


       


for org in CD:
    for con in CD[org]:
        for prots in CD[org][con]:
            clust = []

            for prot in prots:
                if prot.bcd0 == "C":
                    clust.append("C")
                if prot.bcd0 == "D":
                    clust.append("D")
                if prot.bcd1 == "C":
                    clust.append("C")
                if prot.bcd1 == "D":
                    clust.append("D")
                if prot.bcd1 == "CD":
                    clust.append("C")
                    clust.append("D")
            clustNoDups = []
            for let in clust:
                try:
                    clustNoDups.index(let)
                except:
                    clustNoDups.append(let)
            clustNoDups.sort()
            if clustNoDups == ["C", "D"]:
                try:
                    TOMMorgsCD[org]
                    TOMMorgsCD[org].append(prots)
                except:
                    TOMMorgsCD[org] = [prots]
                

    
for d in [TOMMorgsBCD, TOMMorgsCD]:
    for org in d:
        for prots in d[org]:
            consolidate_TOMM(prots)



for org in TOMMorgsBCD:
    eliminate_TOMMBCD(TOMMorgsBCD[org], TOMMorgsCD)

    
for org in TOMMorgsCD:
    eliminate_TOMMCD(TOMMorgsCD[org])

    

    
"""organize by protein type"""
DprotsBCD =[]
CprotsBCD = []
BprotsBCD = []
CDprotsBCD = []
other = []
for org in TOMMorgsBCD:
    if TOMMorgsBCD[org] != []:
        for prots in TOMMorgsBCD[org]:
            check_B(prots)
            for prot in tuple(prots):
                if prot.orgQuery == "Bacillus sp. B14905":
                    continue
                if prot.bcd1 == "D" or (prot.bcd0 == "D" and prot.bcd1 == "?"):
                    DprotsBCD.append(prot)
                if prot.bcd1 == "C" or (prot.bcd0 == "C" and prot.bcd1 == "?"):
                    CprotsBCD.append(prot)
                if prot.bcd1 == "B"or (prot.bcd0 == "B" and prot.bcd1 == "?")or prot.bcd1 == "?BC":
                    BprotsBCD.append(prot)
                if prot.bcd1 == "CD":
                    CDprotsBCD.append(prot)
                if prot.bcd1 != "D" and prot.bcd1 != "C" and prot.bcd1 != "B" and prot.bcd1 != "CD":
                    other.append(prot)


for lis in [DprotsBCD, CprotsBCD, BprotsBCD, CDprotsBCD, other]:
    File_writer(lis[0].bcd1 + "protsBCD", lis)

    
DprotsCD =[]
CprotsCD = []
CDprotsCD = []
other = []
for org in TOMMorgsCD:
    if TOMMorgsCD[org] != []:
        for prots in TOMMorgsCD[org]:
            
            for prot in prots:
                if prot.bcd1 == "D" or (prot.bcd0 == "D" and prot.bcd1 == "?"):
                    DprotsCD.append(prot)
                if prot.bcd1 == "C" or (prot.bcd0 == "C" and prot.bcd1 == "?"):
                    CprotsCD.append(prot)
                if prot.bcd1 == "CD":
                    CDprotsCD.append(prot)
                if prot.bcd1 != "D" and prot.bcd1 != "C" and prot.bcd1 != "B" and prot.bcd1 != "CD":
                    other.append(prot)


    
for lis in [DprotsCD, CprotsCD, CDprotsCD]:
    File_writer(lis[0].bcd1 + "protsCD", lis) #CprotsCD comes out as ?protsCD.txt instead....fix this later

"""For Writing Files based on organism"""

for d in [TOMMorgsBCD, TOMMorgsCD]:
    for org in d:
        if d[org] != []:
            org_file_writer(org, d[org])

def overhang(s, e, leng):
    s = str(int(s))
    e = str(int(e) +20)
    if int(s) < 21:
        s = "1"
    if int(e) > (leng - 20):
        e = str(leng)
    return s, e

f = open("ChalfCD_BCD.txt", "w")
for prot in CDprotsBCD:
    f.write(">" + "gi" + str(prot.gid) + "|" + prot.org + "|" + prot.bcd0 + "|" + prot.bcd1 + "|" + "\n")
    f.write(prot.aa_seq[:len(prot.aa_seq)/2] + "\n")
    f.write("\n")
f.close()


f = open("DhalfCD_BCD.txt", "w")
for prot in CDprotsBCD:
    f.write(">" + "gi" + str(prot.gid) + "|" + prot.org + "|" + prot.bcd0 + "|" + prot.bcd1 + "|" + "\n")
    f.write(prot.aa_seq[-len(prot.aa_seq)/2:] + "\n")
    f.write("\n")
f.close()

f = open("DhalfCD_CD.txt", "w")
for prot in CDprotsCD:
    f.write(">" + "gi" + str(prot.gid) + "|" + prot.org + "|" + prot.bcd0 + "|" + prot.bcd1 + "|" + "\n")
    f.write(prot.aa_seq[-len(prot.aa_seq)/2:] + "\n")
    f.write("\n")
f.close()

f = open("ChalfCD_CD.txt", "w")
for prot in CDprotsCD:
    f.write(">" + "gi" + str(prot.gid) + "|" + prot.org + "|" + prot.bcd0 + "|" + prot.bcd1 + "|" + "\n")
    f.write(prot.aa_seq[:len(prot.aa_seq)/2] + "\n")
    f.write("\n")
f.close()

f = open("mcbC_regions2.1.txt", "w")
g = open("mcbC_regions2.2.txt", "w")
h = open("mcbC_regions1.txt", "w")

for prot in BprotsBCD:    
    gfile = prot.gbfile
    switch = False
    count = 0
    for line in tuple(gfile):
        if line.find("FEATURES") != -1:
            switch = True
        if switch:
            if line.find('/region_name="mcbC-like_oxidoreductase"') != -1:
                count += 1
    switch = False
    count2 = 1
    for line in tuple(gfile):
        if line.find("FEATURES") != -1:
            switch = True
        if switch:
            if line.find('/region_name="mcbC-like_oxidoreductase"') != -1:
                if count > 1:
                    if count2 == 1:
                        reg = gfile[gfile.index(line)-1]
                        s, e = re.split(" +", reg[:-1])[2].split("..")
                        s, e = re.search("[0-9]+", s).group(), re.search("[0-9]+", e).group()
                        s, e = overhang(s, e, prot.length)
                        f.write(">" + "gi" + str(prot.gid) + "|" + prot.org + "|" + prot.bcd0 + "|" + prot.bcd1 + "|" + str(count2) + "\n")
                        f.write(prot.aa_seq[int(s)-1:int(e)] + "\n")
                        f.write("\n")
                        count2 += 1
                        gfile.remove(line) # WARNING!!!! THIS WILL REMOVE THIS LINE FORM THE PROTEIN'S GBFILE ATTRIBUTE!!!! DO NOT SAVE PROTEINS WITH THIS PART OF SCRIPT ON!!!!
                        continue
                    if count2 == 2:
                        reg = gfile[gfile.index(line)-1]
                        s, e = re.split(" +", reg[:-1])[2].split("..")
                        s, e = re.search("[0-9]+", s).group(), re.search("[0-9]+", e).group()
                        s, e = overhang(s, e, prot.length)
                        g.write(">" + "gi" + str(prot.gid) + "|" + prot.org + "|" + prot.bcd0 + "|" + prot.bcd1 + "|" + str(count2) + "\n")
                        g.write(prot.aa_seq[int(s)-1:int(e)] + "\n")
                        g.write("\n")
                        gfile.remove(line) # WARNING!!!! THIS WILL REMOVE THIS LINE FORM THE PROTEIN'S GBFILE ATTRIBUTE!!!! DO NOT SAVE PROTEINS WITH THIS PART OF SCRIPT ON!!!!
                if count == 1:
                    reg = gfile[gfile.index(line)-1]
                    s, e = re.split(" +", reg[:-1])[2].split("..")
                    s, e = re.search("[0-9]+", s).group(), re.search("[0-9]+", e).group()
                    s, e = overhang(s, e, prot.length)
                    h.write(">" + "gi" + str(prot.gid) + "|" + prot.org + "|" + prot.bcd0 + "|" + prot.bcd1 + "|" + str(count2) + "\n")
                    h.write(prot.aa_seq[int(s)-1:int(e)] + "\n")
                    h.write("\n")
                    gfile.remove(line) # WARNING!!!! THIS WILL REMOVE THIS LINE FORM THE PROTEIN'S GBFILE ATTRIBUTE!!!! DO NOT SAVE PROTEINS WITH THIS PART OF SCRIPT ON!!!!                 


f.close()
g.close()
h.close()


def ycao_split(prots, cfile, dfile, CDTorF):
    f = open(dfile, "w")
    if CDTorF:
        g = open(cfile, "w")
    
    for prot in prots:    
        gfile = prot.gbfile
        switch = False
        
        for line in tuple(gfile):
            if line.find("FEATURES") != -1:
                switch = True
            if switch:
                if line.find('/region_name="YcaO"') != -1:
                    reg = gfile[gfile.index(line)-1]
                    s, e = re.split(" +", reg[:-1])[2].split("..")
                    s, e = re.search("[0-9]+", s).group(), re.search("[0-9]+", e).group()
                    s, e = overhang(s, e, prot.length)
                    f.write(">" + "gi" + str(prot.gid) + "|" + prot.org + "|" + prot.bcd0 + "|" + prot.bcd1 + "|" + str(count2) + "\n")
                    f.write(prot.aa_seq[int(s)-1:] + "\n")
                    f.write("\n")
                    if CDTorF:
                        g.write(">" + "gi" + str(prot.gid) + "|" + prot.org + "|" + prot.bcd0 + "|" + prot.bcd1 + "|" + str(count2) + "\n")
                        g.write(prot.aa_seq[:int(s)+20] + "\n")
                        g.write("\n")
    f.close()
    if CDTorF:
        g.close()

for prots, cfile, dfile, CDTorF in [
    [CDprotsBCD, "CD_cycloBCD.txt", "CD_ycaoBCD.txt", True],
    [DprotsBCD, "", "D_ycaoBCD.txt", False],
    [CDprotsCD, "CD_cycloCD.txt", "CD_ycaoCD.txt", True],
    [DprotsCD, "", "D_ycaoCD.txt", False]]: ycao_split(prots, cfile, dfile, CDTorF)
    






        
    

                
"""For Listing tomm proteins by organism"""
##for org in TOMMorgsBCD:
##	print org
##	for prots in TOMMorgsBCD[org]:
##		for prot in prots:
##			try:
##				print prot.gid, prot.bcd0, prot.bcd1
##			except:
##				print prot.gid, prot.bcd0, prot.bcd1
##				break
