# Benchmarking of Methylationcaller for Nanopore Sequencing Data
Benchmarking of the Methylationcaller for Nanopore Sequencing data (Bachelor thesis Marvin Mayer)


## Snakefile
This Shell script is able to run the methylationcalling for the callers deepsignal and tombo. For the input basecalled multi fast5 files are needed which will be unpacked into single fast5 files in the first step. Afterwards the resquiggle of the single fast5 files is done. Both methylationcaller need resquiggled single fast5 files as an input for the methylationcalling.

After methylation values are gathered the output is analysed in two different ways. First the methylation status for each CpG is calculated and summed up into the input format for methrix (Chromosome; Position; number of methylated reads; number of unmethylated reads). 
Second the outputs of the methylationcaller is compared to the true dataset of a bisulfite sequencing. Therefore all the positions which are present in the bisulfite file and in the methylation file is written into a new comparison matrix.

## Snakefile_centre_comparison

To compare the results of the different centres with eachother and with the bisulfite file, this Snakefile is needed. It has the same functionalities as the Snakefile before with the addition of comparing the results for each centre.
