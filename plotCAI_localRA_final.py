import sys,re,os,getopt
import numpy as np
from numpy import mean
from Bio import SeqIO
from scipy import stats
import matplotlib as mpl
mpl.use('Agg') # enables saving graphics
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from collections import OrderedDict
import math

###############
# SUBROUTINES #
###############


# PDCUB scores file looks like this:
# ENST00000641515.2       OR4F5   60      59      0       2618    2.9214349680077714
def readPdcubScores(pdcubScores):
    pdcubScoresDict = {}
    f1 = open(pdcubScores)
    dataf1 = f1.readlines()
    linePattern = re.compile('(\S*)\s(\S*)\s(\S*)\s(\S*)\s(\S*)\s(\S*)\s(\S*)')
    for line in dataf1:
        if linePattern.search(line):
            match = linePattern.search(line)
            shortID = match.group(1)
            pdcubScore = match.group(7)
            pdcubScoresDict[shortID] = pdcubScore
    f1.close()
    return pdcubScoresDict

def validTranscript(cdsLength,cds,stop):
    if cdsLength % 3 == 0 and cdsLength >= (binCount*binSizeNt) and cds.startswith("AUG") and stop in stops:
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
        geneID = record.id # looks like this: ENST00000641515.2|ENSG00000186092.6|OTTHUMG00000001094|OTTHUMT00000003223.1|OR4F5-202|OR4F5|2618|UTR5:1-60|CDS:61-1041|UTR3:1042-2618|
        shortMatch = geneRecordPattern.search(geneID)
        shortID = shortMatch.group(1)
        if cdsPattern.search(description):
            match = cdsPattern.search(description)
            start = int(match.group(1))-1
            end = int(match.group(2))-1
            cds = transcript[start:end+1]
            cdsLength = len(cds)
            stop = transcript[end-2:end+1]
            if validTranscript(cdsLength,cds,stop):
                if shortID in pdcubScoresDict.keys():
                    valid_mRNAs.append((shortID,geneID,str(transcript),str(cds),start,end,cdsLength))
    return valid_mRNAs

def getGC(valid_mRNAs):
    firstSegmentNt = binSizeNt*[0]
    firstSegmentGc = binSizeNt*[0]
    firstSegmentNtSingleRatios = []
    singleCdsGcRatiosSeg1 = {} # single-value total GC ratio for each individual transcript, averaged across first 150 nt positions
    secondSegmentNt = binSizeNt*[0]
    secondSegmentGc = binSizeNt*[0]
    secondSegmentNtSingleRatios = []
    singleCdsGcRatiosSeg2 = {} # single-value total GC ratio for each individual transcript, averaged across second 150 nt positions
    totalNt = binSizeNt*binCount*[0]
    totalGc = binSizeNt*binCount*[0]
    totalNtSingleRatios = []
    singleCdsGcRatiosAll = {} # single-value total GC ratio for each individual transcript, averaged across first 300 nt positions

    for mRNA in valid_mRNAs:
        shortID,geneID,transcript,cds,start,end,cdsLength = mRNA
        firstSegmentNtSingle = 0
        firstSegmentGcSingle = 0
        secondSegmentNtSingle = 0
        secondSegmentGcSingle = 0
        totalNtSingle = 0
        totalGcSingle = 0

        # first segment GC ratio
        for position,nucleotide in enumerate(cds[0:binSizeNt]):
            if nucleotide:
                firstSegmentNt[position] += 1
                firstSegmentNtSingle += 1
            if nucleotide in ["G","C"]:
                firstSegmentGc[position] += 1
                firstSegmentGcSingle += 1
        firstSegmentGcSingleRatio = float(firstSegmentGcSingle)/firstSegmentNtSingle
        singleCdsGcRatiosSeg1[shortID] = firstSegmentGcSingleRatio

        # second segment GC ratio
        for position,nucleotide in enumerate(cds[binSizeNt:binSizeNt*binCount]):
            if nucleotide:
                secondSegmentNt[position] += 1
                secondSegmentNtSingle += 1
            if nucleotide in ["G","C"]:
                secondSegmentGc[position] += 1
                secondSegmentGcSingle += 1
        secondSegmentGcSingleRatio = float(secondSegmentGcSingle)/secondSegmentNtSingle
        singleCdsGcRatiosSeg2[shortID] = secondSegmentGcSingleRatio

        # total GC ratio
        for position,nucleotide in enumerate(cds[0:binSizeNt*binCount]):
            if nucleotide:
                totalNt[position] += 1
                totalNtSingle += 1
            if nucleotide in ["G","C"]:
                totalGc[position] += 1
                totalGcSingle += 1
        totalGcSingleRatio = float(totalGcSingle)/totalNtSingle
        singleCdsGcRatiosAll[shortID] = totalGcSingleRatio

    singleNtGcRatiosSeg1 = [float(i)/float(j) for i,j in zip(firstSegmentGc,firstSegmentNt)] # GC ratio for each nt position, for first 150 nt, averaged across all transcripts
    singleNtGcRatiosSeg2 = [float(i)/float(j) for i,j in zip(secondSegmentGc,secondSegmentNt)] # GC ratio for each nt position, for second 150 nt, averaged across all transcripts
    singleNtGcRatiosAll = [float(i)/float(j) for i,j in zip(totalGc,totalNt)] # GC ratio for each nt position, for first 300 nt, averaged across all transcripts
    return singleNtGcRatiosSeg1,singleCdsGcRatiosSeg1,singleNtGcRatiosSeg2,singleCdsGcRatiosSeg2,singleNtGcRatiosAll,singleCdsGcRatiosAll


def getCodonFrequencies(valid_mRNAs):
    codonBinOne = []
    codonBinTwo = []
    codonBinTotal = []
    codonBinPositional = {}
    for codonPosition in range(0,binSizeNt*binCount/3):
        codonBinPositional[codonPosition] = []

    codonFrequenciesOne = {}
    codonFrequenciesTwo = {}
    codonFrequenciesTotal = {}
    codonFrequenciesPerPosition = {}
    codonFrequenciesPerPositionCumulative = {}
    
    for mRNA in valid_mRNAs:
        shortID,geneID,transcript,cds,start,end,cdsLength = mRNA
        for i in range(0,binSizeNt,3):
            codon = cds[i:i+3]
            codonBinOne.append(str(codon))
        for i in range(binSizeNt,binSizeNt*binCount,3):
            codon = cds[i:i+3]
            codonBinTwo.append(str(codon))
        for i in range(0,binSizeNt*binCount,3): #binSizeNt*binCount,3):
            codon = cds[i:i+3]
            codonBinTotal.append(str(codon))
            codonBinPositional[i/3].append(str(codon))

    for codon in codonMap:
        codonFrequenciesOne[codon] = codonBinOne.count(codon)
        
        codonFrequenciesTwo[codon] = codonBinTwo.count(codon)
        
        codonFrequenciesTotal[codon] = codonBinTotal.count(codon)

        codonFrequenciesPerPosition[codon] = {}
        for codonPosition in range(binSizeNt*binCount/3):
            codonFrequenciesPerPosition[codon][codonPosition] = codonBinPositional[codonPosition].count(codon)

        codonFrequenciesPerPositionCumulative[codon] = {}
        positionalSum = 0
        for codonPosition in range(binSizeNt*binCount/3):
            positionalSum = positionalSum + codonFrequenciesPerPosition[codon][codonPosition]
            codonFrequenciesPerPositionCumulative[codon][codonPosition] = positionalSum

    return codonFrequenciesOne,codonFrequenciesTwo,codonFrequenciesTotal,codonFrequenciesPerPosition,codonFrequenciesPerPositionCumulative

def calculateRelativeAdaptiveness(codonFrequenciesOne,codonFrequenciesTwo,codonFrequenciesTotal,codonFrequenciesPerPosition,codonFrequenciesPerPositionCumulative):
    RA1group = {}
    RA2group = {}
    RAtotalgroup = {}
    cumulativePositionalRAgroup = {}
    for position in range(0,codonSequenceLength):
        cumulativePositionalRAgroup[position] = {}

    # calculate first bin
    for codon in codonFrequenciesOne:
        synonymGroup = residueMap[codonMap[codon]]
        synonymCounts = {}
        for synonym in synonymGroup:
            if synonym in codonFrequenciesOne.keys():
                synonymCounts[synonym] = codonFrequenciesOne[synonym]
        maxSynonymCount = max(synonymCounts.values())
        if maxSynonymCount >> 0:
            RA1 = float(codonFrequenciesOne[codon]) / float(maxSynonymCount)
            RA1group[codon] = RA1

    # calculate second bin
    for codon in codonFrequenciesTwo:
        synonymGroup = residueMap[codonMap[codon]]
        synonymCounts = {}
        for synonym in synonymGroup:
            if synonym in codonFrequenciesTwo.keys():
                synonymCounts[synonym] = codonFrequenciesTwo[synonym]
        maxSynonymCount = max(synonymCounts.values())
        if maxSynonymCount >> 0:
            RA2 = float(codonFrequenciesTwo[codon]) / float(maxSynonymCount)
            RA2group[codon] = RA2

    # calculate total
    for codon in codonFrequenciesTotal:
        synonymGroup = residueMap[codonMap[codon]]
        synonymCounts = {}
        for synonym in synonymGroup:
            if synonym in codonFrequenciesTotal.keys():
                synonymCounts[synonym] = codonFrequenciesTotal[synonym]
        maxSynonymCount = max(synonymCounts.values())
        if maxSynonymCount >> 0:
            RAtotal = float(codonFrequenciesTotal[codon]) / float(maxSynonymCount)
            RAtotalgroup[codon] = RAtotal


    # calculate cumulative positional
    # codonFrequenciesPerPositionCumulative looks like this:
    # {'ACC': {0: 0, 1: 982, 2: 2294, 3: 3313, 4: 4272, 5: 5301 ....
    for codon in codonFrequenciesPerPositionCumulative:
        synonymGroup = residueMap[codonMap[codon]]
        synonymPositionalCumulativeCounts = {}
        for position in codonFrequenciesPerPositionCumulative[codon].keys():
            synonymPositionalCumulativeCounts[position] = {}
            for synonym in synonymGroup:
                if synonym in codonFrequenciesPerPositionCumulative.keys():
                    synonymPositionalCumulativeCounts[position][synonym] = codonFrequenciesPerPositionCumulative[synonym][position]
            positionalMaxSynonymCount = max(synonymPositionalCumulativeCounts[position].values())
            if positionalMaxSynonymCount >> 0:
                cumulativePositionalRA = float(codonFrequenciesPerPositionCumulative[codon][position]) / float(positionalMaxSynonymCount)
                cumulativePositionalRAgroup[position][codon] = cumulativePositionalRA

    return RA1group, RA2group, RAtotalgroup, cumulativePositionalRAgroup

def calculateCAIs(valid_mRNAs,RA1group,RA2group,RAtotalgroup,cumulativePositionalRAgroup):

    firstSegmentCAIs = {}
    secondSegmentCAIs = {}
    totalCAIs = {}
    transcriptCumulativePositionalCodons = {}
    cumulativePositionalCAIs = {}
    cumulativeAveragePositionalCAIs = {}

    for mRNA in valid_mRNAs:
        firstSegmentCodons = []
        secondSegmentCodons = []
        totalCodons = []
        shortID,geneID,transcript,cds,start,end,transcriptLength = mRNA
        transcriptCumulativePositionalCodons[shortID] = {}
        for i in range(0,binSizeNt,3):
            codon = cds[i:i+3]
            firstSegmentCodons.append(str(codon))
        for i in range(binSizeNt,2*binSizeNt,3):
            codon = cds[i:i+3]
            secondSegmentCodons.append(str(codon))
        for i in range(0,2*binSizeNt,3):
            transcriptCumulativePositionalCodons[shortID][i/3] = []
            codon = cds[i:i+3]
            totalCodons.append(str(codon))
            transcriptCumulativePositionalCodons[shortID][i/3].append(str(codon))
            if i>0:
                transcriptCumulativePositionalCodons[shortID][i/3] = transcriptCumulativePositionalCodons[shortID][(i-1)/3] + transcriptCumulativePositionalCodons[shortID][i/3]
        # totalCodons will look like this:
        # ['AUG', 'CCU', 'GAG', 'AUC', 'AGA', ....

        # first partial segment of transcript
        codonRA1values = []
        for codonName in firstSegmentCodons:
            codonRA1values.append(RA1group[codonName])
        firstSegmentCodonsProduct = float(np.prod(codonRA1values))
        firstSegmentCodonLength = float(len(firstSegmentCodons))
        firstSegmentCAI = firstSegmentCodonsProduct**(1/firstSegmentCodonLength)
        firstSegmentCAIs[shortID] = firstSegmentCAI

        # second partial segment of transcript
        codonRA2values = []
        for codonName in secondSegmentCodons:
            codonRA2values.append(RA2group[codonName])
        secondSegmentCodonsProduct = float(np.prod(codonRA2values))
        secondSegmentCodonLength = float(len(secondSegmentCodons))
        secondSegmentCAI = secondSegmentCodonsProduct**(1/secondSegmentCodonLength)
        secondSegmentCAIs[shortID] = secondSegmentCAI

        # full length
        # RAtotalgroup looks like this:
        # {'ACC': 1.0, 'GUC': 0.4918708260688582, 'ACA': 0.7672590136598513, ....
        codonRAtotalvalues = []
        for codonName in totalCodons:
            codonRAtotalvalues.append(RAtotalgroup[codonName])
        totalCodonsProduct = float(np.prod(codonRAtotalvalues))
        totalCodonLength = float(len(totalCodons))
        totalCAI = totalCodonsProduct**(1/totalCodonLength)
        totalCAIs[shortID] = totalCAI

        # cumulative positional
        # transcriptCumulativePositionalCodons looks like this:
        # {'ENST00000494801.5': {0: ['AUG'], 1: ['AUG', 'GAA'], ....
        # cumulativePositionalRAgroup looks like this:
        # {0: {'AUG': 1.0}, 1: {'ACC': 1.0, 'GUC': 0.38607115821347465, ....
        cumulativePositionalCAIs[shortID] = {}
        for position in range(0,codonSequenceLength):
            codonRAcumulativepositionalvalues = []
            for codonName in transcriptCumulativePositionalCodons[shortID][position]:
                codonRAcumulativepositionalvalues.append(cumulativePositionalRAgroup[position][codonName])
            cumulativePositionalCodonsProduct = float(np.prod(codonRAcumulativepositionalvalues))
            positionalCodonLength = float(len(transcriptCumulativePositionalCodons[shortID][position]))
            cumulativePositionalCAI = cumulativePositionalCodonsProduct**(1/positionalCodonLength)
            cumulativePositionalCAIs[shortID][position] = cumulativePositionalCAI
    for position in range(0,codonSequenceLength):
        for shortID in cumulativePositionalCAIs.keys():
            cumulativeAveragePositionalCAIs[position] = np.mean(cumulativePositionalCAIs[shortID][position])


    return firstSegmentCAIs,secondSegmentCAIs,totalCAIs,cumulativeAveragePositionalCAIs

        
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

########
# MAIN #
########

usage = "Usage: python " + sys.argv[0] + " <fastafile> <PDCUB score file> <CDS window size in nt (multiple of 3)>"

if len(sys.argv) != 4:
    print usage
    sys.exit()
    
fastaFile = sys.argv[1]
pdcubScores = sys.argv[2]
windowSizeNt = sys.argv[3]

binSizeNt = int(windowSizeNt)
binCount = 2 # number of bins to examine
codonSequenceLength = binCount*binSizeNt/3 # this is total examined length of sequence, counting by codons

pdcubScoresDict = readPdcubScores(pdcubScores)
valid_mRNAs = processFastaFile(fastaFile)
singleNtGcRatiosSeg1,singleCdsGcRatiosSeg1,singleNtGcRatiosSeg2,singleCdsGcRatiosSeg2,singleNtGcRatiosAll,singleCdsGcRatiosAll = getGC(valid_mRNAs)
# singleNtGcRatiosSeg1 contains individual GC ratio for each nt position, for first 150 nt, averaged across all transcripts
# singleCdsGcRatiosSeg1 contains single-value total GC ratio for each individual transcript, averaged across first 150 nt positions    
# similarly for variables ending in '2' or 'All'

singleCodonGcRatiosAll = []
for position in range (0,len(singleNtGcRatiosAll),3):
    if position > len(singleNtGcRatiosAll)-3:
        break
    tripletSum = float(singleNtGcRatiosAll[position]+singleNtGcRatiosAll[position+1]+singleNtGcRatiosAll[position+2])
    tripletAverage = tripletSum/3
    singleCodonGcRatiosAll.append(tripletAverage)

codonFrequenciesOne,codonFrequenciesTwo,codonFrequenciesTotal,codonFrequenciesPerPosition,codonFrequenciesPerPositionCumulative = getCodonFrequencies(valid_mRNAs)
RA1group,RA2group,RAtotalgroup,cumulativePositionalRAgroup = calculateRelativeAdaptiveness(codonFrequenciesOne,codonFrequenciesTwo,codonFrequenciesTotal,codonFrequenciesPerPosition,codonFrequenciesPerPositionCumulative)
firstSegmentCAIs,secondSegmentCAIs,totalCAIs,cumulativeAveragePositionalCAIs = calculateCAIs(valid_mRNAs,RA1group,RA2group,RAtotalgroup,cumulativePositionalRAgroup)

# plot CAI of second segment against CAI of first segment
plt.figure(dpi=600)

colors = ['midnightblue','lime']
colorMap = LinearSegmentedColormap.from_list('dah', colors, N=100)

for gene in totalCAIs:
    if totalCAIs[gene] < 0.1:
        print gene,totalCAIs[gene]

caiSegOneVals = []
caiSegTwoVals = []
caiTotalVals = []
gcPerCdsSegOneVals = []
gcPerCdsSegTwoVals = []
gcPerCdsAllVals = []
cols = []
IDs = []
for ID in firstSegmentCAIs.keys():
    IDs.append(ID)
for ID in IDs:
    caiSegOneVals.append(float(firstSegmentCAIs[ID]))
    caiSegTwoVals.append(float(secondSegmentCAIs[ID]))
    caiTotalVals.append(float(totalCAIs[ID]))
    gcPerCdsSegOneVals.append(float(singleCdsGcRatiosSeg1[ID]))
    gcPerCdsSegTwoVals.append(float(singleCdsGcRatiosSeg2[ID]))
    gcPerCdsAllVals.append(float(singleCdsGcRatiosAll[ID]))
    cols.append(float(pdcubScoresDict[ID]))


scatterPoints = list(zip(caiSegOneVals,caiSegTwoVals,caiTotalVals,gcPerCdsSegOneVals,gcPerCdsSegTwoVals,gcPerCdsAllVals,cols,IDs))
scatterPoints.sort(key=lambda x:x[6],reverse=False) # True means low values on top of plot
caiSegOneVals,caiSegTwoVals,caiTotalVals,gcPerCdsSegOneVals,gcPerCdsSegTwoVals,gcPerCdsAllVals,cols,IDs = zip(*scatterPoints)

plt.scatter(caiSegOneVals,caiSegTwoVals,c=cols,s=1,cmap=colorMap,edgecolor='none',alpha=0.8)

slope, intercept, r_value, p_value, std_err = stats.linregress(caiSegOneVals,caiSegTwoVals)
r_squared = r_value*r_value
def bestfitGraph(regressionFormula, x_range):
    x = np.array(x_range)
    y = regressionFormula(x)
    plt.plot(x,y,c='blue',label='{:.2}'.format(slope)+'x + '+'{:.2}'.format(intercept)+', $R^2$={:.2}'.format(r_squared),lw=1.2,linestyle='dashed')
def regressionFormula(x):
    return x*slope + intercept
bestfitGraph(regressionFormula, [x/10.0 for x in range(2,11,1)])

plt.colorbar(label='PDCUB score')

#for index,codonName in enumerate(RA2codons):
#    plt.annotate(codonName,(RA1values[index],RA2values[index]),fontsize=5)
#    if codonName in TICs:
#        plt.annotate(codonName,(RA1values[index],RA2values[index]),fontsize=5,color='red')
plt.plot([0.4,1],[0.4,1],label="y = x",c='black',lw=0.3)
plt.xlim([0.4,1.0])
plt.ylim([0.4,1.0])
plt.legend(framealpha=0.5)
plt.xlabel("CAI (0 to " + str(binSizeNt) + " nt)")
plt.ylabel("CAI ( " + str(binSizeNt) + " to " + str(binSizeNt*binCount) + " nt)")
plt.title("position-dependent CAI")
plt.savefig("posDepCAI_localRA_" + str(binSizeNt)  + "bin.pdf")
plt.savefig("posDepCAI_localRA_" + str(binSizeNt)  + "bin.png")
plt.close()
