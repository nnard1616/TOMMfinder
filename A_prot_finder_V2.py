from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re

##orgList = ["Brachyspira pilosicoli 95/1000",
##           "Brachyspira murdochii DSM 12563",
##           "Brachyspira hyodysenteriae WA1",
##           "Borrelia spielmanii A14S",
##           "Staphylococcus lugdunensis HKU09-01",
##           "Actinosynnema mirum DSM 43827",
##           "Bacillus cereus ATCC 14579",
##           "Propionibacterium acnes KPA171202",
##           "Propionibacterium acnes SK137",
##           "Propionibacterium acnes SK187",
##           "Frankia sp. CcI3",
##           "Lysinibacillus sphaericus C3-41",
##           "Bacillus cereus NVH0597-99",
##           u'Thermobispora bispora DSM 43833',
##           u'Rhodococcus erythropolis PR4',
##           u'Streptomyces scabiei 87.22',
##           u'Micromonospora sp. ATCC 39149',
##           u'Rhodococcus erythropolis SK121',
##           u'Streptomyces clavuligerus ATCC 27064',
##           u'Nonomuraea sp. WU8817',
##           u'Frankia symbiont of Datisca glomerata',
##           u'Clavibacter michiganensis subsp. sepedonicus',
##           u'Actinomyces odontolyticus F0309',
##           u'Actinomadura melliaura',
##           u'Salinispora arenicola CNS-205',
##           u'Nonomuraea sp. Bp3714-39',
##           "Corynebacterium urealyticum DSM 7109"]
orgList = ["Streptomyces sp. C"]

orgDNAdict = {}

class DNA_seq:
	def __init__(self, dna1, dna2, dna3, dnaRC1, dnaRC2, dnaRC3, start, end):
		self.dna1 = dna1
		self.dna2 = dna2
		self.dna3 = dna3
		self.dnaRC1 = dnaRC1
		self.dnaRC2 = dnaRC2
		self.dnaRC3 = dnaRC3
		self.start = start
		self.end = end

def get_DNA_sequence(prot):
        Entrez.email = "nnard1616@gmail.com"
        handle = Entrez.efetch(db="nucleotide", id=prot.contig, rettype = "fasta")
        seq_record = SeqIO.read(handle, "fasta")
        contigDNA = str(seq_record.seq)
        start = prot.start
        end = prot.end
        if prot.start.find("<") != -1 or prot.start.find(">") != -1:
                start = prot.start[1:]
        else:
                start = prot.start
        if prot.end.find("<") != -1 or prot.end.find(">") != -1:
                end = prot.end[1:]
        else:
                end = prot.end
        end = int(end) + 15000
        start = int(start) - 15000
        if int(end) > len(contigDNA) - 1:
                end = len(contigDNA) -1
        if int(start) < 0:
                start = 0
        
        DNA = contigDNA[int(start):int(end)]
        return DNA, start, end


def get_A(prot, f=""):
##        haha = False
        DNA, start, end = get_DNA_sequence(prot)
        DNA = Seq(DNA)
        DNA = DNA_seq(DNA, DNA[1:], DNA[2:], DNA.reverse_complement(), DNA.reverse_complement()[1:], DNA.reverse_complement()[2:], start, end)
        hits = []
        seqs = [DNA.dna1,DNA.dna2,DNA.dna3, DNA.dnaRC1, DNA.dnaRC2, DNA.dnaRC3]
        orgDNAdict[prot.org] = seqs
        for seq in seqs:
                tranProt = seq.translate()
                frags = str(tranProt).split("*")
                for frag in frags:

                        frag = frag[-130:]

                        if len(frag) <= 130 and len(frag) > 35 and (frag.find("M") != -1 or frag.find("V") != -1):
                                if len(frag) - frag.find("M") >= 35:
                                        ts = frag[-25:].count("T")
                                        cs = frag[-25:].count("C")
                                        ss = frag[-25:].count("S")
                                        GG = len(re.findall("[GA][A-Z]?[SCT]", frag[-30:]))
                                        cstp = re.findall("[CST]+[^CST]{0,4}[CST]+", frag[-30:])
                                        count = ts+cs+ss
                                        
                                        weight = 12.0*cs+10*ss+ts
                                        percent = (ts + cs + ss) / 25.0

##                                        for clus in cstp:
##                                                clusC = clus.count("C")
##                                                clusS = clus.count("S")
##                                                clusT = clus.count("T")
##                                                if(clusC+clusS+clusT) <= 6:
##                                                        cstp.remove(clus)
##                                                        count -= (clusC+clusS+clusT)
##                                                        cs -=clusC
##                                                        ss -=clusS
##                                                        ts -=clusT
                                        
                                        
                                        try:
                                                clustering = count/float(len(cstp))
                                        except:
                                                clustering = 0
                                        score = clustering*weight 
                                        if percent >= 0.0 and weight > 0 and score > 200: #put score limit to 500
##                                                haha = True
                                                print prot.org, prot.gid
                                                print frag, "from" + str(seqs.index(seq)), " of length ", str(len(frag))
                                                print cs, ss, ts, percent, clustering, weight, score
                                                print cstp
                                                print str(tranProt).index(frag)*3
                                                print
        return DNA, seqs
##                                                f.write(">" + prot.org + " | " + prot.gid + "\n")
##                                                f.write(frag + "\n")
##                                                f.write("\n")
##                                        if haha:
##                                                try:
##                                                        missed.index(prot.org)
##                                                except:
##                                                        missed.append(prot.org)
                                        
missed = []
for org in orgList:
        for prots in TOMMorgsBCD[org]:
                DNA, seqs = get_A(prots[0])

##f= open("actinobacteriaPrecursorPeptides3.txt", "w")


##for org in TOMMorgsBCD:
##        for prots in TOMMorgsBCD[org]:
##                for prot in prots:
##                        for line in prot.gbfile:
##                                if line.find("Actinobacteria") != -1:
##                                        print org, prot.gid
##                                        get_A(prot, f)
##                                        print
##                                        break
##                        break
##
##f.close()
