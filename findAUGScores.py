
#computes PDCUB score for all the AUGS in a transcript. Reports whether the start codon has the highest AUG score
#along with other relevant features of the transcript such as length, gene name, etc in a TSV file.

#Author: Kaavya Subramanian


import sys, re
import math
from Bio import SeqIO
import statistics
import csv
import scipy.stats as st
import numpy as np
import random


def createTSV(filename,data,colNames):

    with open(filename,'wb') as f_output:
        tsv_output = csv.writer(f_output, delimiter = '\t')
        tsv_output.writerow(colNames)
        for row in data:
            tsv_output.writerow(row)


def parseDefline(record):

    match  = re.search('CDS:(\d+)-(\d+)', record.id)
    code = str(record.id).split('|')[0]
    gene = str(record.id).split('|')[5]
    start = int(match.group(1))
    end = int(match.group(2))
    transcript = record.seq.transcribe()
    length = len(transcript)

    return start,end,transcript,code,gene,length




def processTSV(tsvFile):

    colNames = []
    data = []
    with open(tsvFile) as t:
        reader = csv.reader(t, delimiter='\t')
        colNames = next(reader)
        colNames = colNames[1:]

        for row in reader:
            data.append(row)


    for i in range(0,len(data)):
        for j in range(1,len(data[i])):
            data[i][j] = float(data[i][j])


    logVals = {}

    for i in range(0, len(data)):
        binVals = {}
        for j in range(1,len(data[i])):
            binVals[j-1] = data[i][j]

        logVals[data[i][0]] = binVals


    return logVals

def processFastaFile(fastaFile):

    #define all information we want to record about the transcript
    starts, ends, transcripts, transcriptIDS, genes, lengths = [],[],[],[],[],[]


    for record in SeqIO.parse(fastaFile, "fasta"):

        start,end,transcript,transcriptID,gene,length = parseDefline(record)

        if len(transcript[start-1:end]) % 3 != 0: continue
        if transcript[start-1:start+2] != 'AUG': continue
        if transcript[end-3:end] not in ['UAA','UAG','UGA']: continue

        #all at the same index
        starts.append(start -1)
        ends.append(end-1)
        transcripts.append(transcript)
        transcriptIDS.append(transcriptID)
        genes.append(gene)
        lengths.append(length)

    return starts, ends, transcripts, transcriptIDS, genes, lengths

def assignScore(AUGcodons, AUGcodonVals, nullOption,startPos):

    score = 0
    scScore = 0.0
    searchIndex = AUGcodons.index(startPos) #identify start codon position with AUGcodons list
    if AUGcodonVals[searchIndex] == max(AUGcodonVals): #start codon has highest AUG score
        score = 1

    scScore = AUGcodonVals[searchIndex]

    if nullOption == "nullModel":  #random null model option

        highestScore = max(AUGcodonVals)
        randomStart = random.choice(AUGcodonVals)

        if randomStart == highestScore: score = 1
        else: score = 0

    return score, scScore

def returnAUGcodonVals_GC(AUGcodons,logData,transcript, numBins,binSize):
        AUGcodonVals = [] #holds all scores of start codons
        for x in AUGcodons:
            iOffset = x+1
            score = 0.0
            for j in range(0,300):
                if j + iOffset >= len(transcript): break
                nuc = transcript[j+iOffset]
                if nuc == 'G' or nuc == 'C':
                    score += logData['S'][j]
                else: score += logData['W'][j]

            AUGcodonVals.append(score)
        return AUGcodonVals

def returnAUGcodonVals_AUGC(AUGcodons,logData,transcript, numBins,binSize):
        AUGcodonVals = [] #holds all scores of start codons
        for x in AUGcodons:
            iOffset = x+1
            score = 0.0
            for j in range(0,300):
                if j + iOffset >= len(transcript): break
                nuc = transcript[j+iOffset]
                score += logData[nuc][j]


            AUGcodonVals.append(score)
        return AUGcodonVals


def returnAUGcodonVals(AUGcodons,logData,transcript,numBins,binSize):

    AUGcodonVals = [] #holds all scores of start codons
    for x in AUGcodons:
        iOffset = x+3
        score = 0.0
        for j in range(0, numBins):
            for k in range(0,binSize,3):
                if iOffset+k+3 > len(transcript): break
                triplet = transcript[iOffset+k:iOffset+k+3]

                if triplet == 'UAG' or triplet == 'UAA' or triplet == 'UGA':
                    score += 0
                else:
                    score += logData[triplet][j]


            iOffset = iOffset+binSize
        AUGcodonVals.append(score)

    return AUGcodonVals


def buildScoreFile(logData,fastaFile,numBins,binSize,modelType):

    tsvData = []

    starts, ends, transcripts, transcriptIDS, genes, lengths = processFastaFile(fastaFile)


    n = 0   #ensures all information can be accessed with same index

    #compute log scores for all start codons within a transcript
    for transcript in transcripts:
        sys.stdout.write("Transcript %d/%d    \r" %(n+1,len(transcripts)))
        sys.stdout.flush()

        AUGcodons = []    #holds all indices of start codons in transcript
        for match in re.finditer('AUG',str(transcript)):
            index = match.start()
            AUGcodons.append(index)


        if modelType == "gc":
            AUGcodonVals = returnAUGcodonVals_GC(AUGcodons,logData,transcript,numBins,binSize)
        elif modelType == "nucleotide":
            AUGcodonVals = returnAUGcodonVals_AUGC(AUGcodons,logData,transcript,numBins,binSize)
        else:
            AUGcodonVals = returnAUGcodonVals(AUGcodons,logData,transcript,numBins,binSize)



        #identify highest start codon
        s, scScore = assignScore(AUGcodons,AUGcodonVals,nullOption,starts[n])

        utr5Len = starts[n] - 1

        tsvData.append([transcriptIDS[n],genes[n],starts[n],utr5Len,s,lengths[n],scScore])
        n += 1
        #if n >= 10000: break

    return tsvData

def main():

    usage = "Usage: " + sys.argv[0] + "<log TSV>" + "<FASTA file>" + "<null option>" + "<bin size>" + "<bin #> + <modelType>"
    if len(sys.argv) != 7:
        print(usage)
        sys.exit()

    logFile = sys.argv[1]   #file of log scores
    fastaFile = sys.argv[2] #file with gene transcripts
    nullOption = sys.argv[3]   #null model option
    binSize = int(sys.argv[4])   #number of nucleotides in a bin
    numBins = int(sys.argv[5])   #number of bins used to compute score
    modelType = sys.argv[6]

    logData = processTSV(logFile)

    tsvData = buildScoreFile(logData,fastaFile,numBins,binSize,modelType)

    fileName = "highestAUGscores_{}_binSize{}_numBins{}.tsv".format(nullOption,binSize,numBins)

    createTSV(fileName,tsvData,["TranscriptID","Gene","Start Codon","5UTR","Start Codon Highest Score","Length", "Score"])


if __name__ == '__main__':
    main()
