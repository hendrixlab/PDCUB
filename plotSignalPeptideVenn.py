import sys,re,os,getopt
from scipy import stats
from scipy.stats import binom_test
import numpy as np
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg') # enables saving graphics
import matplotlib.pyplot as plt
from matplotlib_venn import venn2,venn3
from collections import OrderedDict

###############
# SUBROUTINES #
###############

# retrieve list of gene IDs from qval set
# file1 format looks like this:
# geneID                  protein     length    log score            p-val                      q-val
# ENST00000215832.11      MAPK1       5881      11.08212786285333    2.1972759969529324e-23     7.50809208158817e-20
def readQValueFile(qValueFile):
        qValueDict = {}
        f1 = open(qValueFile)
        dataf1 = f1.readlines()
        linePattern = re.compile('(\S*)') #\s(\S*)\s(\S*)\s(\S*)\s(\S*)\s(\S*)')
        for line in dataf1:
                if linePattern.search(line):
                        match = linePattern.search(line)
                        geneID = match.group(1)
                        qValueDict[geneID] = True
        f1.close()
        return qValueDict

# retrieve gene IDs and SignalP scores for all valid proteins from gencode human mRNA transcriptome set
# file2 format looks like this:
# geneID              SignalP classificaiton     signal peptide score     "other" score
# ENST00000641515.2   OTHER                      0.000304                 0.999696
def readSignalPOutputFile(signalPOutputFile):
        signalpScoreList = []
        f2 = open(signalPOutputFile)
        dataf2 = f2.readlines()
        linePattern = re.compile('(\S*)\s(\S*)\s(\S*)\s(\S*)')
        for line in dataf2:
                if linePattern.search(line):
                        match = linePattern.search(line)
                        geneID = match.group(1)
                        signalpScore = float(match.group(3))
                        if signalpScore > 0.5:
                                signalpScoreList.append( ( geneID, float(signalpScore) ) )
        f2.close()
        return signalpScoreList

# retrieve gene IDs and Predisi scores for all valid proteins from gencode human mRNA transcriptome set
# file2 format looks like this:
# geneID              Predisi classificaiton     signal peptide score     "other" score
# ENST00000641515.2   OTHER                      0.000304                 0.999696
def readPredisiOutputFile(predisiOutputFile):
        predisiScoreList = []
        f2 = open(predisiOutputFile)
        dataf2 = f2.readlines()
        linePattern = re.compile('(\S*)\s(\S*)\s(\S*)\s(\S*)')
        for line in dataf2:
                if linePattern.search(line):
                        match = linePattern.search(line)
                        geneID = match.group(1)
                        predisiScore = float(match.group(2))
                        if predisiScore > 0.5:
                                predisiScoreList.append( ( geneID, float(predisiScore) ) )
        f2.close()
        return predisiScoreList

####################
# GLOBAL VARIABLES #
####################

stops = ['UAA','UGA','UAG']

TICs = ['UAC','AAC','UUC','UAU','AAG','GAG','AAU','GAU','AUC','GAC','GUG']

codonMap = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
            "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
            "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
            "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
            "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
            "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
            "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
            "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
            "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
            "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
            "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
            "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
            "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
            "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
            "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
            "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

residueMap = {"F":("UUU","UUC"), "L":("UUA","UUG","CUU","CUC","CUA","CUG"), "S":("UCU","UCC","UCA","UCG","AGU","AGC"),
              "Y":("UAU","UAC"), "*":("UAA","UGA","UAG"), "C":("UGU","UGC"), "W":("UGG",), "P":("CCU","CCC","CCA","CCG"),
              "H":("CAU","CAC"), "Q":("CAA","CAG"), "R":("CGU","CGC","CGA","CGG","AGA","AGG"), "I":("AUU","AUC","AUA"),
              "M":("AUG",), "T":("ACU","ACC","ACA","ACG"), "N":("AAU","AAC"), "K":("AAA","AAG"), "V":("GUU","GUC","GUA","GUG"),
              "A":("GCU","GCC","GCA","GCG"), "D":("GAU","GAC"), "E":("GAA","GAG"), "G":("GGU","GGC","GGA","GGG")}

##########
## MAIN ##
##########

usage = "Usage: python " + sys.argv[0] + " <q-value file> <signalP output file> <prediSi output file>"

if len(sys.argv) != 4:
        print usage
        sys.exit()

qValueFile = sys.argv[1]
signalPOutputFile = sys.argv[2]
predisiOutputFile = sys.argv[3]

print('Reading q-value file.')
qValueDict = readQValueFile(qValueFile) # retrieve list of gene IDs from qval set
print('Reading signalP file.')
signalpScoreList = readSignalPOutputFile(signalPOutputFile) # retrieve gene IDs and SignalP scores for all valid proteins from gencode human mRNA transcriptome set
print('Reading prediSi file.')
predisiScoreList = readPredisiOutputFile(predisiOutputFile)

signalpSet = [entry[0] for entry in signalpScoreList] # this is SignalP geneIDs for ONLY the intersection of the gencode and qval sets
predisiSet = [entry[0] for entry in predisiScoreList] # this is PrediSi geneIDs for ONLY the intersection of the gencode and qval sets

qvalSet = []
for key in qValueDict:
        qvalSet.append(key)

print('Plotting venn diagrams.')

plt.figure()
venn3([set(signalpSet), set(predisiSet), set(qvalSet)], set_labels = ('SignalP (all GENCODE)','PrediSi (all GENCODE)','significant PDCUB'))
plt.savefig("vennThreeway.pdf")
plt.savefig("vennThreeway.png")
plt.close()

print('Venn plotted.')
