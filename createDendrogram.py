#Takes in a tsv file and plots into a heatmap and dendrogram.
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
import scipy
import pylab
import scipy.cluster.hierarchy as sch

def zDiff(z):
    return float(z[2])-float(z[1])

def parseargs():
    parser = argparse.ArgumentParser()

    parser.add_argument("input_file",help="input TSV data", type=str)
    parser.add_argument("output_file",help="Output PDF File",type=str)
    parser.add_argument("colormap",help="Color scheme",type=str)
    parser.add_argument("colorTitle",help="Title of Colorbar",type=str)
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



def createCorrMatrix(data):
    #creating distance matrix based on sum of squares
    rawdata = []
    for codon in data:
        rawdata.append(codon)


    dMatrix = scipy.zeros([len(rawdata),len(rawdata)])

    for i in range(len(rawdata)):
        for j in range(len(rawdata)):
            s = 0.0
            for k in range(0,len(rawdata[j])):
                sq = (rawdata[i][k] * rawdata[j][k])
                s += sq

            dMatrix[i,j] = s

    return dMatrix

def heatmap(output_file,heatmapType, data,row_labels, col_labels,cmap='',cbarlabel=''):

   #Normalizes the colorbar to zero as the middle threshold.
   norm = MidpointNormalize(midpoint=0)


   fig = pylab.figure(figsize=(12,8.5))


   [ax1_x, ax1_y, ax1_w, ax1_h] = [0.05,0.22,0.2,0.6]
   [axm_x, axm_y, axm_w, axm_h] = [0.285,0.9,0.4,0.6]
   axm_y = ax1_y
   axm_h = ax1_h

   [axcb_x, axcb_y,axcb_w,axcb_h] = [0.55, 0.22, 0.2, 0.6]


   dMatrix = createCorrMatrix(data)
   fig = pylab.figure()
   ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=True)
   Y = sch.linkage(dMatrix, method='average')
   Z = sch.dendrogram(Y,orientation='left', color_threshold=None)
   ax1.axis('off')

   axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h],frame_on=True)
   newdata = np.array(data)
   idx1 = Z['leaves']
   newdata = newdata[idx1]
   row_labels = np.array(row_labels)
   row_labels = row_labels[idx1]

   im = axm.matshow(newdata, aspect='auto', origin='lower',cmap=cmap, norm=norm)


   ax = fig.add_axes([axcb_x,axcb_y, axcb_w, axcb_h], frame_on=False)
   fig.colorbar(im, cmap=cmap, ax=ax, norm=norm, orientation='vertical',label=cbarlabel)
   ax.set_xticks([])
   ax.set_yticks([])


   #Set and label the x and y axes.

   start, end = axm.get_xlim()
   axm.xaxis.set_ticks(np.arange(start,end,10))
   axm.set_yticks(np.arange(len(row_labels)))
   axm.set_xticklabels(np.arange(0,3000,600))
   axm.set_xlabel("Position after Start Codon (nt)", fontsize = 10)
   axm.set_yticklabels(row_labels,fontsize = 6)
   axm.tick_params(axis='x', which='major', labelsize=8)
   axm.tick_params(axis='y', which='major', labelsize=4)

   axm.tick_params(top=False,bottom=True, labeltop=False,labelbottom=True)


   fig.tight_layout()
   pylab.savefig(output_file)



#parses TSV file for data and labels to create heatmap
def processTSV(tsvfile):

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



    rowNames = [row[0] for row in data]
    data = [row[1:] for row in data]


    return rowNames,colNames,data

def main():
    args = parseargs()
    rowNames,colNames,data = processTSV(args.input_file)
    heatmap(args.output_file,heatmapType,data,rowNames,colNames,cmap=args.colormap,cbarlabel = args.colorTitle)


if __name__ == '__main__':
    main()
