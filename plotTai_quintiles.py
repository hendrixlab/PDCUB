import sys,re,os,getopt
import numpy as np
from numpy import mean
from Bio import SeqIO
from scipy import stats
from scipy.stats import binom_test
from array import *
import matplotlib as mpl
mpl.use('Agg') # enables saving graphics
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib_venn import venn2,venn3
from matplotlib.colors import LinearSegmentedColormap
from collections import OrderedDict
import math

###############
# SUBROUTINES #
###############

def validTranscript(cdsLength,cds,initiator,stop):
    if cdsLength % 3 == 0 and cdsLength >= 3*windowSize and initiator == "AUG" and stop in stops:
        return True
    else:
        return False

def processFastaFile(fastaFile):
    cdsPattern = re.compile('CDS:(\d+)\-(\d+)')
    geneRecordPattern = re.compile('(ENST\d*\.\d*)\S*')
    sequences = SeqIO.parse(fastaFile,"fasta")
    valid_mRNAs = []
    for record in sequences:
        description = record.description
        transcript = record.seq
        transcript = transcript.upper()
        transcript = transcript.transcribe()
        geneID = record.id
        # geneID looks like: ENST00000641515.2|ENSG00000186092.6|OTTHUMG00000001094|OTTHUMT00000003223.1|OR4F5-202|OR4F5|2618|UTR5:1-60|CDS:61-1041|UTR3:1042-2618|
        shortMatch = geneRecordPattern.search(geneID)
        shortID = shortMatch.group(1)
        if cdsPattern.search(description):
            match = cdsPattern.search(description)
            start = int(match.group(1))-1
            end = int(match.group(2))-1
            cds = transcript[start+3:end-2]
            cdsLength = len(cds)
            initiator = transcript[start:start+3]
            stop = transcript[end-2:end+1]
            if validTranscript(cdsLength,cds,initiator,stop):
                valid_mRNAs.append((shortID,geneID,str(transcript),str(cds),start,end,cdsLength))
    return valid_mRNAs

def extractTaiData(valid_mRNAs):
    avgTaiScores = [ 0.0 ] * totalScores
    totalTaiScores = [ 0 ] * totalScores

    for mRNA in valid_mRNAs:
        shortID,geneID,transcript,cds,start,end,cdsLength = mRNA
        taiScores = [ 1.0 ] * totalScores

        # w is in units of codons, w defines a window
        for w in range(0,totalScores):

            # wStart is the actual window start position in nt
            wStart = w*3
            if len(cds) >= wStart + 3*windowSize: 
                # i is the position of a codon in units of nucleotides
                for i in range(wStart,wStart+3*windowSize,3):
                    codon = cds[i:i+3]
                    taiScores[w] *= tAI[codon]

                # take the n-th root for this window w 
                taiScores[w] = taiScores[w]**(1.0/windowSize)

                # add this particular transcript's tAI scores to the total tAI sum for all transcripts
                avgTaiScores[w] += taiScores[w]
                totalTaiScores[w] += 1

    xValues = []
    yValues = []
    for w in range(totalScores):
        #print(w, avgTaiScores[w]/totalTaiScores[w])
        xValues.append(w)
        yValues.append(avgTaiScores[w]/totalTaiScores[w])
    return xValues,yValues
         
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

aaType = {}
aaValue = {}
aaType["R"] = "+"
aaValue["R"] = 'Blues'
aaType["H"] = "+"
aaValue["H"] = 'Blues'
aaType["K"] = "+"
aaValue["K"] = 'Blues'
aaType["D"] = "-"
aaValue["D"] = 'Greens'
aaType["E"] = "-"
aaValue["E"] = 'Greens'
aaType["S"] = "P"
aaValue["S"] = 'Purples'
aaType["T"] = "P"
aaValue["T"] = 'Purples'
aaType["N"] = "P"
aaValue["N"] = 'Purples'
aaType["Q"] =  "P"
aaValue["Q"] = 'Purples'
aaType["Y"] =  "P"
aaValue["Y"] = 'Purples'
aaType["C"] = "N"
aaValue["C"] = 'Oranges'
aaType["G"] =  "N"
aaValue["G"] = 'Oranges'
aaType["P"] =  "N"
aaValue["P"] = 'Oranges'
aaType["A"] =  "N"
aaValue["A"] = 'Oranges'
aaType["I"] =  "N"
aaValue["I"] = 'Oranges'
aaType["L"] =  "N"
aaValue["L"] = 'Oranges'
aaType["M"] =  "N"
aaValue["M"] = 'Oranges'
aaType["F"] =  "N"
aaValue["F"] = 'Oranges'
aaType["W"] =  "N"
aaValue["W"] = 'Oranges'
aaType["V"] =  "N"
aaValue["V"] = 'Oranges'

# tAI scores for codons in humans (scores retrieved from Tuller 2010 translational ramp paper)
tAI = {}
tAI["UUU"] = 0.161002
tAI["UUC"] = 0.366748
tAI["UUA"] = 0.213936
tAI["UUG"] = 0.282396
tAI["UCU"] = 0.336186
tAI["UCC"] = 0.242054
tAI["UCA"] = 0.152845
tAI["UCG"] = 0.171149
tAI["UAU"] = 0.218399
tAI["UAC"] = 0.449878
tAI["UAA"] = 0.061128
tAI["UAG"] = 0.050122
tAI["UGU"] = 0.402506
tAI["UGC"] = 0.91687
tAI["UGA"] = 0.091687
tAI["UGG"] = 0.304401
tAI["CUU"] = 0.366748
tAI["CUC"] = 0.264059
tAI["CUA"] = 0.091724
tAI["CUG"] = 0.334963
tAI["CCU"] = 0.305623
tAI["CCC"] = 0.220049
tAI["CCA"] = 0.213967
tAI["CCG"] = 0.190709
tAI["CAU"] = 0.147586
tAI["CAC"] = 0.336186
tAI["CAA"] = 0.336186
tAI["CAG"] = 0.749389
tAI["CGU"] = 0.213936
tAI["CGC"] = 0.154034
tAI["CGA"] = 0.183395
tAI["CGG"] = 0.211491
tAI["AUU"] = 0.535208
tAI["AUC"] = 0.552567
tAI["AUA"] = 0.152855
tAI["AUG"] = 0.611247
tAI["ACU"] = 0.305623
tAI["ACC"] = 0.220049
tAI["ACA"] = 0.183405
tAI["ACG"] = 0.242054
tAI["AAU"] = 0.459902
tAI["AAC"] = 1
tAI["AAA"] = 0.519563
tAI["AAG"] = 0.685819
tAI["AGU"] = 0.107335
tAI["AGC"] = 0.244499
tAI["AGA"] = 0.183374
tAI["AGG"] = 0.211491
tAI["GUU"] = 0.336186
tAI["GUC"] = 0.242054
tAI["GUA"] = 0.152845
tAI["GUG"] = 0.537897
tAI["GCU"] = 0.886308
tAI["GCC"] = 0.638142
tAI["GCA"] = 0.27515
tAI["GCG"] = 0.240831
tAI["GAU"] = 0.254921
tAI["GAC"] = 0.580685
tAI["GAA"] = 0.397311
tAI["GAG"] = 0.52445
tAI["GGU"] = 0.201253
tAI["GGC"] = 0.458435
tAI["GGA"] = 0.275061
tAI["GGG"] = 0.301956

########
# MAIN #
########

usage = "Usage: python " + sys.argv[0] + " <fastafile1> <fastafile2> <fastafile3> <fastafile4> <fastafile5> <CDS window size in codons>"

if len(sys.argv) != 7:
    print usage
    sys.exit()
    
fastaFile1 = sys.argv[1]
fastaFile2 = sys.argv[2]
fastaFile3 = sys.argv[3]
fastaFile4 = sys.argv[4]
fastaFile5 = sys.argv[5]
windowSize = int(sys.argv[6])

totalScores = 400

valid_mRNAs1 = processFastaFile(fastaFile1)
valid_mRNAs2 = processFastaFile(fastaFile2)
valid_mRNAs3 = processFastaFile(fastaFile3)
valid_mRNAs4 = processFastaFile(fastaFile4)
valid_mRNAs5 = processFastaFile(fastaFile5)

xValues1,yValues1 = extractTaiData(valid_mRNAs1)
xValues2,yValues2 = extractTaiData(valid_mRNAs2)
xValues3,yValues3 = extractTaiData(valid_mRNAs3)
xValues4,yValues4 = extractTaiData(valid_mRNAs4)
xValues5,yValues5 = extractTaiData(valid_mRNAs5)

plt.figure(figsize=[6.4,3.7],dpi=600)
ax = plt.axes()
label1='1 (lowest)'
label2='2 (20%-40%)'
label3='3 (40%-60%)'
label4='4 (60%-80%)'
label5='5 (highest)'
plt.plot(xValues1,yValues1,c="#fa9fb5",label=label1)
plt.plot(xValues2,yValues2,c="#f768a1",label=label2)
plt.plot(xValues3,yValues3,c="#dd3497",label=label3)
plt.plot(xValues4,yValues4,c="#ae017e",label=label4)
plt.plot(xValues5,yValues5,c="#7a0177",label=label5)
handles, labels = ax.get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ordered_label_list = [label5,label4,label3,label2,label1]
ordered_label_values = [by_label[k] for k in ordered_label_list]
lg = ax.legend(ordered_label_values, ordered_label_list, framealpha=0.9, fontsize=6, title='PDCUB quintiles', ncol=1, loc='lower right')
lg.get_title().set_fontsize(6)
plt.title('tAI of human transcriptome by PDCUB score')
plt.xlabel('Distance from start codon (codons)')
plt.ylabel('Local tAI')
plt.savefig('taiWindow' + str(windowSize) + '.quintiles.' + str(totalScores) + '.pdf')
plt.savefig('taiWindow' + str(windowSize) + '.quintiles.' + str(totalScores) + '.png')
plt.close()
