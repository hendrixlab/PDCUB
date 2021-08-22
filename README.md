# PDCUB
A collection of tools for computing position-dependent codon usage bias (PDCUB) scores.


### computeCodonStats.py
Performs a position-based computational survey for each non-stop codon given a FASTA file. 

Required input:
* GENCODE gene FASTA
* Path to the result directory (must already exist)
* bin size in nucleotides, must be a multiple of 3
* CDS cutoff length in nucleotides (the maximum CDS position to read to)

Example: `python computeCodonStats.py gene.fa result_dir 60 600`

Script output:
* TSV files where the rows are codons, the columns are the bins, and the values are
  * codon raw counts per bin
  * codon z-scores per bin
  * codon p-values per bin
  * codon z-squared values per bin
  * codon log likelihood weights per bin
 
### computeCodonStats_GC.py
Performs a similar computation as computeCodonStats.py, but computes for two categories:
* G/C nucleotides, or S nucleotides
* A/U nucleotides, or W nucleotides

Inputs are the same as computeCodonStats.py, with the omission of a bin size.

Example: `python computeCodonStats_GC.py gene.fa result_dir 600`

### computeCodonStats_AUGC.py
Performs a similar computation as computeCodonStats.py, but computes for each nucleotide instead of each codon.

Inputs are the same as computeCodonStats.py, with the omission of a bin size.

Example: `python computeCodonStats_AUGC.py gene.fa result_dir 600`

### findAUGScores.py
Finds the highest PDCUB scoring AUG and identifies whether it is the start codon for each transcript in an input FASTA file.

Required input:
* log-likelihood TSV that was computed from computeCodonStats.py, computeCodonStats_GC.py, or computeCodonStats_AUGC.py
* gene FASTA file 
* option whether or not to choose a random AUG start codon to compare to the highest scoring AUG
  - enter "nullmodel" to generate a null model
  - enter "realmodel" to generate a PDCUB model
* bin size that the log-likelihood TSV was calculated
* the number of bins that the log-likelihood TSV was calculated with
* the type of model to generate, depending on what computing script was used
 - enter "gc" if the log likelihood TSV came from computeCodonStats_GC.py
 - enter "augc" if the log likelihood TSV came from computeCodonStats_AUGC.py
 - enter "codon" if the log likelihood TSV came from computeCodonStats_GC.py

example: `python findAUGScores.py log_gc.tsv gene.fa realmodel 3 100 gc`

Script Output:
* TSV file with the following columns
 * Transcript ID
 * Gene ID
 * Start Codon position 
 * 5' UTR position
 * Whether or not the start codon was the highest scoring AUG
  - 1 if the start codon was the highest scoring AUG
  - 0 if the start codon was not the highest scoring AUG
 * Transcript Length in nucleotides
 * the score of the start codon

### createHeatMap.py
Creates a PDF of a heatmap that is based on one of the computational TSV outputs from one of the computing scripts. In the paper, this script is used to visualize the z-scores of each non-stop codon across the CDS.

Required input:
* input TSV file, where the rows are codons, the columns are bin ranges, and the values are numeric in nature
* name of the output file
* name for the colormap used for the heatmap [colormap options] (https://matplotlib.org/stable/tutorials/colors/colormaps.html)
* The title of the color bar
* Organizational scheme for the heatmap
 - 0: the codons are organized by amino acid and chemical group
 - 1: the codons are organized by chemical group, amino acid, and the first column in increasing order
 - 2: the codons are organized by the difference between the first and second column in increasing order

Example: `python createHeatMap.py z-scores.tsv zScoreHeatMap.pdf plasma z-scores 0`

### createDendrogram.py
Creates a PDF of a heatmap where the codons are organized according to a dendrogram.

The input is the same as createHeatMap.py, but there is no option for the organization scheme of the heatmap.

Example: `python createDendrogram.py z-scores.tsv zScoreDendrogram.py Reds z-scores`

### plotTai_quintiles.py
This python script reads in five evenly sized fasta files, each comprising one quintile of the GENCODE protein-coding human transcriptome sorted by PDCUB score. The output is a combined plot of the average local tAI scores for each quintile. These averages are calculated using a sliding window, the size which is entered from the command line when executing the script.
