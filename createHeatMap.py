#Takes in a tsv file and plots into a heatmap.

#Author: Kaavya Subramanian


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import csv
from matplotlib.colors import Normalize
import sys
import argparse
from operator import itemgetter

def zDiff(z):
    return float(z[2])-float(z[1])

def parseargs():
    parser = argparse.ArgumentParser()

    parser.add_argument("input_file",help="input TSV data", type=str)
    parser.add_argument("output_file",help="Output PDF File",type=str)
    parser.add_argument("colormap",help="Color scheme",type=str)
    parser.add_argument("colorTitle",help="Title of Colorbar",type=str)
    parser.add_argument("heatmapType",help="Organization scheme",type=int,choices=range(0,3))
    args = parser.parse_args()
    return args




#class created by Sahil Kothya from semicolonworld.com
#Class that sets the middle of the colorbar to zero

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

#Function adding a line to the heatmap
#Input given as (x1,y1),(x2,y2), positions relative to the image axes
#lw = line width
def addLine(ax,x1,x2,y1,y2,lw):
    line = plt.Line2D([x1,x2],[y1,y2], transform=ax.transAxes, color='black',lw = lw)
    line.set_clip_on(False)
    ax.add_line(line)



#Function that adds chemical group labels
def createCodonLabelsChem(codonNamesChem,codonNamesAlpha,lw,ax):

    sortedChemList = []

    for codon in codonNamesChem:
        for ids in codon[1]:
            sortedChemList.append(codonNamesAlpha[ids])

    createCodonLabels(sortedChemList,lw,ax)

    topPos = 1
    botPos = 0
    midPos = 0
    counter = 0
    incs = [['+',0.16393],['-',.065574],['polar',.26229],['nonpolar',.50820]]

    for l in incs:
        counter += 1
        botPos = topPos - l[1]
        midPos = ((topPos + botPos) / 2) -.01

        if counter == len(incs):
            botPos = 0

        #Turned off for GOldwater 2020
        #addLine(ax,-.25,-.25,botPos,topPos,lw)
        #addLine(ax,-.25,-.2,topPos,topPos,lw)
        #addLine(ax,-.25,-.2,botPos,botPos,lw)
        #ax.text(-.32,midPos,l[0],fontsize=10,rotation=90,fontweight = 'light',transform=ax.transAxes,clip_on=False)
        addLine(ax,-.1,-.1,botPos,topPos,lw)
        addLine(ax,-.1,-.05,topPos,topPos,lw)
        addLine(ax,-.1,-.05,botPos,botPos,lw)
        ax.text(-.15,midPos,l[0],fontsize=10,rotation=90,fontweight = 'light',transform=ax.transAxes,clip_on=False)

        topPos = botPos



def createCodonLabels(codonNames,lw,ax):
   topPos = 1
   botPos = 0
   midPos = 0
   inc = .01639344262
   counter = 0
   for codon in codonNames:
       counter += 1
       labelRange = codon[1] * inc
       botPos = topPos - labelRange
       midPos = ((topPos + botPos) / 2) -.01

       if counter == len(codonNames):
           botPos = 0

       addLine(ax,-.15,-.15,botPos,topPos,lw)
       addLine(ax,-.15,-.1,topPos,topPos,lw)
       addLine(ax,-.15,-.1,botPos,botPos,lw)
       ax.text(-.2,midPos,codon[0],fontsize=7,fontweight = 'light',transform=ax.transAxes,clip_on=False)

       topPos = botPos




def heatmap(heatmapType, data,row_labels, col_labels, ax=None,cbar_kw={},cbarlabel='', **kwargs):

   #Normalizes the colorbar to zero as the middle threshold.
   norm = MidpointNormalize(midpoint=0)
   im = ax.imshow(data,norm=norm, **kwargs)

   #Create a colorbar.
   cbar = ax.figure.colorbar(im, ax=ax,  **cbar_kw)
   cbar.ax.set_ylabel(cbarlabel, rotation = -90, va="bottom", fontsize = 15)

   #Set and label the x and y axes.
   ax.set_xticks(np.arange(len(col_labels)))
   ax.set_yticks(np.arange(len(row_labels)))
   ax.set_xticklabels(col_labels)
   ax.set_yticklabels(row_labels)

   ax.tick_params(axis='x', which='major', labelsize=8)
   ax.tick_params(axis='y',which='major',labelsize=5)
   ax.tick_params(axis='both', which='minor', labelsize=4)
   ax.set_xticklabels(np.arange(0,3000,300),fontsize = 10)
   ax.set_xlabel("Distance from Start Codon (nt)", fontsize = 12)
   ax.tick_params(top=False,bottom=True, labeltop=False,labelbottom=True)
   #plt.setp(ax.get_xticklabels(),rotation= -90, ha = "right")
   # ax.set_xlabel("Distance from Start Codon (nt)", fontsize = 15)
   # ax.set_ylabel("Codon Z-score Profile", fontsize=15)
   # ax.tick_params(top=False,bottom=False, labeltop=False,labelbottom=False)

   #Label y axis with amino acid labels
   if heatmapType != 2:
      aaNames = [['A',4],['C',2],['E',2],['D',2],['G',4],['F',2],['I',3],['H',2],['K',2],['M',1],['L',6],['N',2],['Q',2],['P',4],['S',6],['R',6],['T',4],['W',1],['V',4],['Y',2]]

      aaNamesChem = [['+',[15,7,8]],['-',[3,2]],['polar',[14,16,11,12,19]],['neutral',[1,4,13,0,6,10,9,5,17,18]]]
      createCodonLabelsChem(aaNamesChem,aaNames,0.3,ax)

   return im,cbar


#parses TSV file for data and labels to create heatmap
def processTSV(tsvfile,heatmapType):
    aaOrdering = ['R','H','K','D','E','S','T','N','Q','Y','C','G','P','A','I','L','M','F','W','V']
    aaDict = {'A':4,'C':2,'E':2,'D':2,'G':4,'F':2,'I':3,'H':2,'K':2,'M':1,'L':6,'N':2,'Q':2,'P':4,'S':6,'R':6,'T':4,'W':1,'V':4,'Y':2}

    residueMap = {"F":("UUU","UUC"), "L":("UUA","UUG","CUU","CUC","CUA","CUG"), "S":("UCU","UCC","UCA","UCG","AGU","AGC"),
                  "Y":("UAU","UAC"), "*":("UAA","UGA","UAG"), "C":("UGU","UGC"), "W":("UGG",), "P":("CCU","CCC","CCA","CCG"),
                  "H":("CAU","CAC"), "Q":("CAA","CAG"), "R":("CGU","CGC","CGA","CGG","AGA","AGG"), "I":("AUU","AUC","AUA"),
                  "M":("AUG",), "T":("ACU","ACC","ACA","ACG"), "N":("AAU","AAC"), "K":("AAA","AAG"), "V":("GUU","GUC","GUA","GUG"),
                  "A":("GCU","GCC","GCA","GCG"), "D":("GAU","GAC"), "E":("GAA","GAG"), "G":("GGU","GGC","GGA","GGG")}
    rowNames = []
    colNames = []
    data = []


    with open(tsvfile) as t:
        reader = csv.reader(t, delimiter = '\t')
        colNames = next(reader)
        colNames = colNames[1:]     #first column is codons, not z-scores


        for row in reader:
           data.append(row)

    #Turn all numbers into float data type
    for i in range(0,len(data)):
        for j in range(1,len(data[i])):
            data[i][j] = float(data[i][j])

    if  heatmapType == 0 or heatmapType == 1:
        reOrderData = []
        for a in aaOrdering:
            for x in range(0,len(data)):
                if data[x][0] in residueMap[a]:
                    y = aaDict[a]
                    for r in range(x, x+y):
                        reOrderData.append(data[r])
                    break

        data = reOrderData


        if heatmapType == 1: #Sorts each amino acid in terms of increasing first bin
            pos = 0
            for n in aaOrdering:
                r = aaDict[n]
                s = data[pos:pos+r]
                s.sort(key=itemgetter(1))
                data[pos:pos+r] = s
                pos = pos + r
    elif heatmapType == 2: #Sorts all rows in terms of difference between first and second bin
        data.sort(key = zDiff)

    rowNames = [row[0] for row in data]
    data = [row[1:] for row in data]


    return rowNames,colNames,data

def main():
    args = parseargs()

    #types of heatmaps:
    #0 = heatmap organized by amino acid and chemical group
    #1 = heatmap organized by chem group, amino acid, and first
    #2 = heatmap organized by difference across first two bins
    heatmapType = args.heatmapType

    rowNames,colNames,data = processTSV(args.input_file,heatmapType)
    fig,ax = plt.subplots()
    im,cbar = heatmap(heatmapType,data,rowNames,colNames,ax=ax,cmap=args.colormap,cbarlabel = args.colorTitle)
    fig.tight_layout()
    plt.savefig(args.output_file)


if __name__ == '__main__':
    main()
