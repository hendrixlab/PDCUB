#Accepts a Gene FASTA file and calculates bin-based GC percentages, z-scores, p-values, and log likelihoods.
#Outputs each result in a TSV file.

#Author: Kaavya Subramanian




import sys, re
import math
from Bio import SeqIO
import statistics
import csv
import scipy.stats as st
#######################################################################################
#Class that represents a bin that the transcript is divided into. Holds information
#like the number of codons in that bin, and the zscore. Each codon has a list of bin
#objects that encompass the entire transcript.
#######################################################################################

class binObject:

    def __init__(self, nameofBin,nucCounts, zscore,percent):
            self.nameofBin = nameofBin      #useful when creating the heat map
            self.nucCounts = nucCounts
            self.zscore = zscore
            self.percent = percent

    def calcZscore(self, mean, sd):
        if sd == 0:     #when there are no codons in the transcript at that position
            self.zscore = 0
        else:
            self.zscore = (self.percent - mean) / sd


    def calcpercent(self,totalCounts):
        #self.percent = self.nucCounts/float(totalCounts)
        self.percent = self.nucCounts/float(totalCounts)

#######################################################################
#Creates a list of Bin objects that every codon will use.
#######################################################################

def createBinList(maxLength):
    numberofBins = maxLength

    binList = []

    for x in range(0, numberofBins):
        nameofBin = "{}".format(x + 1)
        newBin = binObject(nameofBin,0,0,0) #initialize zscore and codon count values as 0
        binList.append(newBin)

    return binList



####################################################################
#Creates a dictionary where the key is the codon and is associated
#with a list of bins.
###################################################################

def createCodonDictionary(maxLength):
    wsValues = {'W': 0, 'S': 0}
    for nuc in wsValues.keys():
        wsValues[nuc] = createBinList(maxLength)

    return wsValues


#Calculates the percents for all codons in the dictionary.

def calcAllPercents(wsValues,allBinValues):
    for key in wsValues.keys():
        n = 0
        for x in wsValues[key]:
            x.calcpercent(allBinValues[n])
            n += 1
#Calculates the z-score for all codons in the dictionary

def calczscores(wsValues):
    for key in wsValues.keys():
        print(key)
        newData = []      #each codon has a unique mean and sd.
        for x in wsValues[key]:
            newData.append(x.percent)

        mean,sd = 0.0,0.0
        if len(newData) > 1:    #there must be at least two values to compute the mean and sd.
            mean = statistics.mean(newData)
            sd = statistics.stdev(newData)
            print(mean, sd)
        for x in wsValues[key]:
            x.calcZscore(mean,sd)

#Sequence Classes and Functions. Returns the CDS of a given transcript based on headerline information.

def parseDefline(record):
    match  = re.search('CDS:(\d+)-(\d+)', record.id)
    start = int(match.group(1))
    end = int(match.group(2))
    rnaSeq = record.seq.transcribe()
    cds = rnaSeq[start-1:end]
    cdsRange = len(cds)
    return cds



#Checks to see if a CDS has any anamolies. Only valid CDS transcripts are returned
def checkGivenCodingSequence(geneTranscript):
    startCodon = "AUG"

    transcriptStartCodon = geneTranscript[0:3]
    transcriptStopCodon = geneTranscript[len(geneTranscript)-3:len(geneTranscript)]


    if (transcriptStartCodon == startCodon) and (transcriptStopCodon in stops) and (len(geneTranscript) % 3 == 0):
        return True
    else:
        return False

#Parses transcripts in FASTA file and extracts CDS. Filters out valid CDSs.
def returnValidSequences(fastaFile, allWSValues):
    print("in valid sequences")
    longestCDS = 0
    validSequences = []

    for record in SeqIO.parse(fastaFile, "fasta"):
        currentSequence = parseDefline(record)
        check = checkGivenCodingSequence(currentSequence)

        if check:


            rnaSeq = record.seq.transcribe()
            allWSValues['W'] += rnaSeq.count('A') + rnaSeq.count('U')
            allWSValues['S'] += rnaSeq.count('G') + rnaSeq.count('C')

            currentSequence = currentSequence[3:]  #exclude start codon
            validSequences.append(currentSequence)
            if len(currentSequence) >= longestCDS:
                longestCDS = len(currentSequence)  #Keeps track of the longest CDS to create the necessary amount of bins.
            else:
                 continue

    return validSequences, longestCDS


def countCodonsintoBins(validSequences,maxLength,wsValues,allBinValues,allWSValues):
    num = 0     #Keeps track of the nth transcript read.
    for x in validSequences:
        num += 1
        sys.stdout.write("Transcript %d/%d    \r" %(num,len(validSequences)))
        sys.stdout.flush()

        for y in range(0,len(x)):
            if y < maxLength:
                nuc = x[y]
                binIndex = y
                if nuc == 'G' or nuc == 'C': wsValues['S'][binIndex].nucCounts += 1
                else: wsValues['W'][binIndex].nucCounts += 1

                allBinValues[binIndex] += 1


def createTSV(filename,data,getBins):
    with open(filename,'wb') as f_output:
        tsv_output = csv.writer(f_output, delimiter = '\t')
        tsv_output.writerow(getBins)
        for row in data:
            tsv_output.writerow(row)


def main():

    usage = "Usage: " + sys.argv[0] + " <Gencode FASTA>" + "<results Directory>" + "<CDS cutoff length (nt)>"
    if len(sys.argv) != 4:
    print(usage)
    sys.exit()

    fastaFile = sys.argv[1]
    resultsDir = sys.argv[2]
    maxLength = int(sys.argv[3])

    stops = ['UAA','UGA','UAG']
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






    allWSValues = {'W':0,'S':0}             #Dictionary that keeps track of the total number of W/S nucleotides
                                            # W = A/U, S = G/C


    validSequences, longestCdsLengths  =  returnValidSequences(fastaFile, allWSValues)

    #Create Data Structures to hold codon values.
    wsValues = createCodonDictionary(maxLength)


    allBinValues = {}
    for b in range(0,maxLength):
        allBinValues[b] = 0

    print("Counting codons into their bins.")
    countCodonsintoBins(validSequences,maxLength,wsValues,allBinValues,allWSValues)




    print("Calcuating percentages.")
    calcAllPercents(wsValues,allBinValues)


    print("Calculating z score vals.")
    calczscores(wsValues)

    #Creates list of Codon names and bin ranges as labels.
    nucNames = ['W','S']

    getBins = []
    getBins.append('Nucleotide')
    for x in wsValues['W']:
        getBins.append(x.nameofBin)


    #Creates arrays that can be read into TSVs.
    zScoreData = []
    for x in nucNames:
        newList = []
        newList.append(x)
        for y in wsValues[x]:
            val = 0 + y.zscore
            newList.append(val)
        zScoreData.append(newList)

    zSquareData = []
    for x in codonNames:
        newList = []
        newList.append(x)
        for y in wsValues[x]:
            val = y.zscore ** 2
            newList.append(val)
        zSquareData.append(newList)


    print("Calculating p-values.")
    #Uses Zscore data to do it since it is essentially a single function change.


    pvalData = []
    for codonRow in zScoreData:
        newList = []
        newList.append(codonRow[0])
        z_counter = 1
        for z_counter in range(1,len(codonRow)):
            p = st.norm.sf(codonRow[z_counter])    #Finds upper tail probability.
            newList.append(p)
        pvalData.append(newList)

    nucCountData = []
    for x in nucNames:
        newList = []
        newList.append(x)
        for y in wsValues[x]:
            val = 0 + y.nucCounts
            newList.append(val)
        nucCountData.append(newList)


    totalNucs = 0
    for x in allWSValues:
        totalNucs += allWSValues[x]

    logLikelihoodData = []
    for x in nucNames:
        newList = []
        newList.append(x)
        for y in wsValues[x]:
            z = allWSValues[x]/float(totalNucs)
            weight = float(y.percent)/(z)
            val = math.log(weight,2)
            newList.append(val)
        logLikelihoodData.append(newList)


    percentData = []
    for x in nucNames:
        newList = []
        newList.append(x)
        for y in wsValues[x]:
            val = 0 + y.percent
            newList.append(val)
        percentData.append(newList)

    createTSV('{}/GCCounts.tsv'.format(resultsDir),nucCountData,getBins)
    createTSV('{}/zscores.tsv'.format(resultsDir),zScoreData,getBins)
    createTSV('{}/pvalues.tsv'.format(resultsDir),pvalData,getBins)
    createTSV('{}/zsquaredvals.tsv'.format(resultsDir), zSquareData,getBins)
    createTSV('{}/percentdata.tsv'.format(resultsDir), percentData, getBins)
    createTSV('{}/logLikelihoodData.tsv'.format(resultsDir),logLikelihoodData, getBins)


if __name__ == 'main':
    main()
