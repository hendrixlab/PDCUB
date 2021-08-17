import os, sys
import re, argparse, itertools
import pysam
from Bio import SeqIO
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter

WINDOW = 19
POLY = 4
MAX_LEN = 1500
MIN_LEN = -150
#at least one read in range described by MAX/MIN
MIN_READS = 1

def main():
    args = get_args()
   
    sequences = read_fasta_pc(args.fpc) 

    #5th quintile is highest
    quints = []
    quints.append(read_fasta_pc(args.q1))
    quints.append(read_fasta_pc(args.q2))
    quints.append(read_fasta_pc(args.q3))
    quints.append(read_fasta_pc(args.q4))
    quints.append(read_fasta_pc(args.q5))

    size = args.size
    
    pc_x = [x for x in range(MIN_LEN//3, MAX_LEN//3+1, size//3)]
    profiles = []
    profiles_windowed = []
    for i in range(len(quints)):
        profiles.append(list(zip(*sorted(get_pc_counts(args.bpc, sequences, quints[i]).items()))))
        profiles[i][1] = np.array(savgol_filter(to_codons(size,profiles[i][1]), 19, 4))
    
    
    colors = ['#c6dbef','#9ecae1','#6baed6','#3182bd','#08519c']
    
    #golden ratio
    fig, ax = plt.subplots(figsize=[6,3.7])
    
    labels = ['1 (lowest)', '2 (20-40%)', '3 (40-60%)', '4 (60-80%)', '5 (highest)']
    for i in range(len(quints)):
        ax.plot(pc_x,100*profiles[i][1],color=colors[i],label=labels[i])

    ax.set_xlim(0,400)
    plt.legend(title='PDCUB quintiles')
    plt.xlabel('Position from CDS (codons)')
    plt.ylabel('Normalized Ribo-seq Profile (%)')
    fig.savefig(args.o + '.pdf')

#sum reads within codon
def to_codons(size, y):
    avg = [0]*(len(y)//size+1)
    for i in range(0,len(y)-size+1,size):
        avg[i//3] = np.sum(y[i:i+size])
    
    return avg


def get_args():
    parser = argparse.ArgumentParser(description="Plot normalized ribosome profile given quintiles (Protein Coding only)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bpc",type=str,required=True,help="BAM file for protein coding Ribo-seq")
    parser.add_argument("--fpc",type=str,required=True,help="Fasta file for all protein coding transcripts")
    parser.add_argument("--q1",type=str,required=True,help="First quintile fasta")
    parser.add_argument("--q2",type=str,required=True,help="Second quintile fasta")
    parser.add_argument("--q3",type=str,required=True,help="Third quintile fasta")
    parser.add_argument("--q4",type=str,required=True,help="Fourth quintile fasta")
    parser.add_argument("--q5",type=str,required=True,help="Fifth quintile fasta")
    parser.add_argument("--o",type=str,required=True,help="Output file name")
    parser.add_argument("--size",type=str,required=False,default=3)
    parser.add_argument("--inset",type=bool,required=False,default=True)

    return parser.parse_args()


def get_pc_counts(bam_file, sequences, quint):
    f = pysam.AlignmentFile(bam_file)
    total_counts = {k:0.0 for k in range(MIN_LEN,MAX_LEN+1)}
    total_transcripts = 0
    quint_ids = [get_CDS_info(seq.id)[1] for seq in quint]
    #list to dict
    quint_ids = dict(itertools.zip_longest(*[iter(quint_ids)] * 2, fillvalue=''))
    
    for seqObj in sequences:
        CDS_start, transcript_ID = get_CDS_info(seqObj.id)
        reads = 0
        temp_counts = {k:0.0 for k in range(MIN_LEN,MAX_LEN+1)}

        if transcript_ID not in quint_ids:
            continue

        if CDS_start is not None:
            for read in f.fetch(seqObj.id):
                positions = read.get_reference_positions()
                relpos = positions[0] - CDS_start + 15
               
                #reads only increments if within range of MIN/MAX
                if relpos in temp_counts:
                    temp_counts[relpos] += 1
                    reads+=1
            if reads >= MIN_READS:
                total_transcripts+=1
                #normalize by number of reads
                for relpos in temp_counts:
                    temp_counts[relpos] /= reads

                #add temp_counts to total_counts
                for relpos in temp_counts:
                    if relpos in total_counts:
                        total_counts[relpos] += temp_counts[relpos]
                    

    #normalize by number of transcripts
    for relpos in total_counts:
        total_counts[relpos] /= total_transcripts

    f.close()

    return total_counts


def get_CDS_info(defline):
    #cds at pos 8, transcript id at pos 0
    terms = defline.split('|')
    try:
        CDS_start = int(terms[8].split(':')[1].split('-')[0])-1
        transcript_ID = terms[0]
        return (int(terms[8].split(':')[1].split('-')[0])-1, terms[0])
    except:
        return (None, None)


def read_fasta_pc(fasta_file):
    fasta = SeqIO.parse(fasta_file, "fasta")
    pc_sequences = []

    for record in fasta:
        pc_sequences.append(record)

    return pc_sequences


main()

