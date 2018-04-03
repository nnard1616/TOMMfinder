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

date = ''
for item in time.gmtime()[0:3]:
  if len(str(item)) ==1:
    date += '0'
    date += str(item)
  else:
    date += str(item)
    
directory = "./%s/" % date  #directory in which dictionaries will be saved (named after the date).  Dictionaries may not be necessary if we start using the computational cluster.


try:
  os.chdir(directory) #moves pointer to that new directory.
except:
  os.mkdir(directory) #creates that directory if it doesn't already exist.
  os.chdir(directory)# moves there after creating it

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
        
"""PREPARATION OF QUERY SEQUENCES"""
BCDorgs = open("../BCDorgs.txt", "r").readlines() #use [:-2] to get rid of "\r\n" at the end of each org name
CDorgs = open("../CDorgs.txt", "r").readlines() #use [:-2] to get rid of "\r\n" at the end of each org name
ALLorgs = BCDorgs + CDorgs

queriesTEST = [["Bacillus amyloliquefaciens FZB42", "154685202"], ["Borrelia spielmanii A14S", "224984462"]]
queriesD = [["Bacillus amyloliquefaciens FZB42", "154685202"],
            ["Borrelia spielmanii A14S", "224984462"],
            ["Clavibacter michiganensis subsp. michiganensis NCPPB 382", "148273958"],
            ["Clostridium botulinum A3 str. Loch Maree", "170761358"],
            ["Herpetosiphon aurantiacus ATCC 23779", "159897446"],
            ["Propionibacterium acnes KPA171202", "50842348"],
            ["Pseudomonas putida KT2440", "26990931"],
            ["Pyrococcus furiosus DSM 3638", "18976374"],
            ["Streptococcus iniae", "23305056"],
            ["Staphylococcus aureus RF122", "82751105"],
            ["Brachyspira hyodysenteriae WA1", "225621255"],
            ["Brevibacterium linens BL2", "260906288"],
            ["Streptococcus gallolyticus UCN34", "288906342"],
            ["Bacillus cereus ATCC 4342", "229154924"],
            ["Bacillus thuringiensis serovar pondicheriensis BGSC 4BA1", "228930592"],
            ["Clostridium botulinum A str. Hall", "153935261"],
            ["Lactobacillus crispatus MV-3A-US", "262046317"],
            ["Listeria monocytogenes Clip81459","226223714"],
            ["Burkholderia mallei GB8 horse 4***","67642384"],
            ["Acidithiobacillus ferrooxidans ATCC 53993","198283860"],
            ["Legionella pneumophila str. Lens***","54292916"],
            ["Escherichia coli","215272915"],
            ["Bacillus pseudofirmus OF4", "288561754"],
            ["Rhodobacterales bacterium HTCC2083", "254461578"],
            ["Rhodobacter sphaeroides ATCC 17029", "126464461"],
            ["Propionibacterium acnes SK187", "289425699"],
            ["Paracoccus denitrificans PD1222", "119386901"],
            ["Clostridium cellulovorans 743B", "302876433"],
            ["Frankia sp. CcI3", "86742875"],
            ["Frankia sp. CcI3", "86742908"],
            ["Lysinibacillus fusiformis ZC1", "299538160"],
            ["Microcystis aeruginosa NIES-298","158934371"],
            ["Bacillus thuringiensis str. Al Hakam", "118476226"],
            ["Prevotella amnii CRIS 21A-A","307565755"],
            ["Bacillus cereus Rock3-44", "229083583"],
            ["Bacillus atrophaeus 1942", "311069755"]
            ]
queriesC = [["Bacillus amyloliquefaciens FZB42", "154685201"],
            ["Borrelia spielmanii A14S", "224984463"],
            ["Clavibacter michiganensis subsp. michiganensis NCPPB 382", "148273959"],
            ["Clostridium botulinum A3 str. Loch Maree", "170760340"],
            ["Herpetosiphon aurantiacus ATCC 23779", "159897445"],
            ["Propionibacterium acnes KPA171202", "50842349"],
            ["Pyrococcus furiosus DSM 3638", "18976375"],
            ["Streptococcus iniae", "23305055"],
            ["Staphylococcus aureus RF122", "82751106"],
            ["Brachyspira hyodysenteriae WA1", "225621254"],
            ["Brevibacterium linens BL2", "260906289"],
            ["Streptococcus gallolyticus UCN34", "288906343"],
            ["Bacillus thuringiensis serovar pondicheriensis BGSC 4BA1", "228930591"],
            ["Bacillus cereus ATCC 4342", "229154924"],
            ["Clostridium botulinum A str. Hall", "153937765"],
            ["Lactobacillus crispatus MV-3A-US", "262046316"],
            ["Listeria monocytogenes Clip81459","226223713"],
            ["Legionella pneumophila str. Lens***","54292916"],
            ["Acidithiobacillus ferrooxidans ATCC 53993","198283859"],
            ["Legionella pneumophila str. Lens***","54292915"],
            ["Propionibacterium acnes SK137","295129710"],
            ["Sulfolobus acidocaldarius DSM 639","70606350"],
            ["Escherichia coli","215272917"],
            ["Pseudomonas putida KT2440","MKTLRVYNYEILNFDSDPMVFSSSGFTRITEPKIIRALQLIEESQSKYIQHHALEKILRKVHLQPSSAMNFLKSLSIIGEPTDPPFFQHTLIYHDLEISDETKNLLESKHSGSLEIRQYSEYTPHTTRKPTLFVFACTKLSPHSLKANYLNLLVSNPHCGATVGFISGNHFHLTEMHIPAIGNPCAFCTLDRITHYEKLRASQHHWCKIWSFCCSSKLDLPKIHVDELQKSLIIGTIVSFISKLTKAPKSKTTQDQVLLSRTINLETGVTTEDPSVHWPLCQCRGVK"],
            ["Bacillus pseudofirmus OF4", "288561753"],
            ["Rhodobacterales bacterium HTCC2083","254460879"],
            ["Rhodobacter sphaeroides ATCC 17029", "126464462"],
            ["Propionibacterium acnes SK187", "289425694"],
            ["Paracoccus denitrificans PD1222", "119386902"],
            ["Clostridium cellulovorans 743B", "302876429"],
            ["Frankia sp. CcI3", "86742875"],
            ["Frankia sp. CcI3", "86742911"],
            ["Lysinibacillus fusiformis ZC1", "299538161"], 
            ["Microcystis aeruginosa NIES-298","158934371"],
            ["Bacillus thuringiensis str. Al Hakam","118476225"],
            ["Actinomyces odontolyticus F0309", "293189640"],
            ["Prevotella amnii CRIS 21A-A","307565739"],
            ["Bacillus cereus Rock3-44", "229083577"],
            ["Bacillus atrophaeus 1942", "311069754"]
            ]
            #Strept Iniae does not have cxxc motifs in c protein, I cannot find c or b protein for pseudomonas putida. Bacillus cereus ATCC 4342, Burkholderia mallei GB8 horse 4,  is a fusion CD. *** indicates that it is a cd only tomm , no b. Mycobacterium tuberculosis H37R is SUPER WEIRD!!!  The c and d are no where near each other!!!!!!!
queriesB = [["Bacillus amyloliquefaciens FZB42", "154685203"],
            ["Borrelia spielmanii A14S", "224984464"],
            ["Clavibacter michiganensis subsp. michiganensis NCPPB 382", "148273957"],
            ["Clostridium botulinum A3 str. Loch Maree", "170759238"],
            ["Herpetosiphon aurantiacus ATCC 23779", "159897448"],
            ["Propionibacterium acnes KPA171202", "50842350"],
            ["Pyrococcus furiosus DSM 3638", "18976373"],
            ["Streptococcus iniae", "23305054"],
            ["Staphylococcus aureus RF122", "82751109"],
            ["Brachyspira hyodysenteriae WA1", "225621253"],
            ["Brevibacterium linens BL2", "260906287"],
            ["Streptococcus gallolyticus UCN34", "288906341"],
            ["Bacillus cereus ATCC 4342", "229154923"],
            ["Bacillus thuringiensis serovar pondicheriensis BGSC 4BA1", "228930593"],
            ["Clostridium botulinum A str. Hall", "153935433"],
            ["Lactobacillus crispatus MV-3A-US", "262046313"],
            ["Listeria monocytogenes Clip81459","226223712"],
            ["Escherichia coli","215272916"],
            ["Pseudomonas putida KT2440","MNISHDEYLKFLSPEVMDETLAFHAKGNYTIHKATQHLSALHQIPPKDLEKLTGNELKLNTIVAPSKNKEPTLTPLASPHPLRRNSSCERFEARPLHFDIVQALLAPLLTKSPTTYKRGYPSGGALYPIEVFCINLNNKIEQWPTESDALHLLPSSRSLEAHTPSIDIHQLTEAIIPESIDIGSPSLALIYSIYLPKALFKYRYRGYRLSLMEAGSMYMITDLRCKELQLNSRPWSGYTDHQVTKNLNLNPSLFLPVCIQLIG"],
            ["Bacillus sp. B14905","126652796"],
            ["Bradyrhizobium japonicum USDA 110", "27379650"],
            ["Lactobacillus gasseri JV-V03","300362803"],
            ["Lysinibacillus sphaericus C3-41","169829450"],
            ["Methylobacterium nodulans ORS 2060","220924544"],
            ["Methylobacterium sp. 4-46","170742111"],
            ["Nonomuraea sp. WU8817","242129435"],
            ["Roseovarius sp. HTCC2601","114765686"],
            ["Streptomyces albus J1074","239981746"],
            ["Bacillus pseudofirmus OF4", "288561755"],
            ["Rhodobacterales bacterium HTCC2083", "254462302"],
            ["Rhodobacter sphaeroides ATCC 17029", "126464460"],
            ["Propionibacterium acnes SK187", "289425702"],
            ["Paracoccus denitrificans PD1222", "119386900"],
            ["Clostridium cellulovorans 743B", "302876430"],
            ["Frankia sp. CcI3", "86742876"],
            ["Frankia sp. CcI3", "86742907"],
            ["Lysinibacillus fusiformis ZC1", "299538163"],
            ["Bacillus licheniformis ATCC 14580", "163119323"],
            ["Microcystis aeruginosa NIES-298", "158934376"],
            ["Bacillus thuringiensis str. Al Hakam", "118476227"],
            ["Prevotella amnii CRIS 21A-A","307565761"],
            ["Bacillus cereus Rock3-44", "229083579"],
            ["Bacillus atrophaeus 1942","311069752"],
            ["Streptomyces griseus subsp. griseus NBRC 13350", "182438204"]
            ]
#can't find good results for Streptomyces hygroscopicus subsp. jinggangensis and Streptomyces hygroscopicus ATCC 53653



"""Helper Functions for my Main function"""

def get_aa_seq(acc): #acc can either be a gid or accession number of some kind
	Entrez.email = "nnard1616@gmail.com"
	try:
		handle = Entrez.efetch(db="protein", rettype = "fasta", id = acc)
	except:
		return "wwwwww", "wwwwwww", "wwwwww"
	try:
		seq_record = SeqIO.read(handle, "fasta")
	except:
		aa_seq, gid, acc = get_aa_seq(acc)
		return aa_seq, gid, acc
	handle.close()
	aa_seq = str(seq_record.seq) #gives amino acid sequence as string
	m = seq_record.id
	m = m.replace("|", "\n")
	m = m.split()
	gid = m[1]
	acc = m[3]
	return aa_seq, gid, acc

def add_to_dict(orgDict, org, prot): #super complicated.... T_T
    try:
        dictList = orgDict[org] #dictList is the list of proteins that return give the organism key --- dictList sounds like dickless...LOL! ---checking to see if org is already in a dictionary key.  If so, dictList is list of proteins that are returned from the key org query to the 
        gidList = dict() 
        for inst in dictList:
            gidList[inst.gid] = inst.bits
        try:
            bits = gidList[prot.gid]#pass means that there is a gid already in the list
            if prot.bits > bits: #if prot that we are trying to add has a larger bit score than the prot that is already there, then replace the old prot with the new one
                for inst in dictList: #now find the location of the old protein to be replaced so that it can be replaced with the new one with a better bit score.
                    if inst.gid == prot.gid:
                        dictList.remove(inst)
                        dictList.append(prot)
        except:
            dictList.append(prot)
    except:
        orgDict[org] = [prot]
    return orgDict

def get_DNA_sequence(prot):
    Entrez.email = "nnard1616@gmail.com"
    handle = Entrez.efetch(db="nucleotide", id=prot.contig, rettype = "fasta")
    seq_record = SeqIO.read(handle, "fasta")
    contigDNA = str(seq_record.seq)
    try:
        DNA = contigDNA[int(prot.start)-1:int(prot.end)]
    except:
        try:
            start = prot.start[1:]
            end = prot.end[1:]
            DNA = contigDNA[int(start)-1:int(end)]
        except:
            return
    return DNA


"""MAIN FUNCTION"""
def __main__(n, orgDict, queries, bcd0):
        for query in queries[n:]:
                if n == len(queries): #if we have reached the end of the query list of proteins, end the recursive function
                        return orgDict
                else:
                        #Be careful with the aa_seq variable, it could be difficult with defining it as it gets redefined often.
                        print n
                        try:
                                result_handle_xml = NCBIWWW.qblast("blastp", "nr", query[1], hitlist_size = 800)
                                n +=1
                        except: #this exception will act as insurance for faulty internet connections
                                TIme = time.localtime()
                                print "something went wrong for " + query[0] + ", at " + str(TIme[3]) + ":" + str(TIme[4]) + " on " + str(TIme[0]) + "/" + str(TIme[1]) + "/" + str(TIme[2]) + "."
                                orgDict = __main__(n, orgDict, queries, bcd0)
                                return orgDict #the dictionary of organsims to proteins!!!
                        blast_record = NCBIXML.read(result_handle_xml)
                        print "finished BLAST of," + str(query[0])
                # # Right here is where I need to insert code to save xml files of my search results
                        for alignment in blast_record.alignments:
                                for hsp in alignment.hsps:
                                        if hsp.bits > 70:
##                                                pdb.set_trace()
                                                line = alignment.title
                                                switch = True
                                                switch2 = True #this switch will, along with tempOrgsList, be used to prevent equivalent, degenerate GID numbers from being added.
                                                tempOrgsList = []
                                                gid = alignment.hit_id
                                                gid = gid.split("gi|")[1][0:gid.split("gi|")[1].find("|")]
##                                                print gid
                                                aa_seq2, gid, acc = get_aa_seq(gid)
                                                while switch:
                                                        org = line[line.find("gi|")+3:line.find("]")+1]
                                                        orgName = org[org.find("[")+1:org.find("]")]
                                                        try:
                                                            tempOrgsList.index(orgName)
                                                            switch2 = False #there already is an entry for this org
                                                        except:
                                                            tempOrgsList.append(orgName)
                                                            gid = org[0:org.find("|")]

                                                        if (line.find("gi|") != -1) and switch2:
                                                            
                                                            prot = Protein(bcd0, gid, len(aa_seq2), hsp.bits, query[0] )
                                                            prot.aa_seq = aa_seq2
                                                            prot.org = orgName
##                                                            print prot.gid, orgName
                                                            orgDict = add_to_dict(orgDict, orgName, prot)
																
                                                        if line.find("gi|") == -1:
                                                                switch = False
                                                        line = line.replace("gi|" + org, "")
                                                        switch2 = True


                        print len(orgDict) #should give number of organisms
                                                
                        result_handle_xml.close()
        return orgDict

orgDictD = __main__(0, dict(), queriesD, "D")
orgDictC = __main__(0, dict(), queriesC, "C")
orgDictB = __main__(0, dict(), queriesB, "B")

dicts = [orgDictD, orgDictC, orgDictB]

d = {}
d[0] = "D"
d[1] = "C"
d[2] = "B"

for dic in dicts:
    f= open("orgDict" + d[dicts.index(dic)] + ".txt", "w")
    pickle.dump(dic, f)
    f.close()
