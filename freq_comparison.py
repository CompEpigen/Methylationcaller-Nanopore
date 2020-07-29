# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:22:06 2020

@author: mmayer

"""

import sys
import csv
import argparse
import os

def proper_round(num, dec=0):
    num = str(num)[:str(num).index('.')+dec+2]
    if num[-1]>='5':
        return float(num[:-2-(not dec)]+str(int(num[-2-(not dec)])+1))
    return float(num[:-1])


def make_key(c, s):
    return str(c) + ":" + str(s)

def to_int(string):
    return int(float(string))

class MethylationStats:
    def __init__(self, chromosome, position, num_reads, num_methylated):
        self.chr = chromosome
        self.pos = to_int(position)
        self.num_reads = to_int(num_reads)
        self.num_methylated_reads = to_int(num_methylated)

    def accumulate(self, num_reads, num_methylated):
        self.num_reads += to_int(num_reads)
        self.num_methylated_reads += to_int(num_methylated)

    def methylation_frequency(self):
        return float(self.num_methylated_reads) / self.num_reads

    def chromosome(self):
        return self.chr

    def position(self):
        return self.pos

    def methylated_reads(self):
        return self.num_methylated_reads

    def reads(self):
        return self.num_reads

    def unmethylated(self):
        unmethylated = self.num_reads - self.num_methylated_reads
        return unmethylated

"""
deepsignal format:
    chrom | pos | strand | pos_in_strand | prob_0_sum | prob_1_sum | count_modified | count_unmotifed
    coverage | modification_frequency | k_mer
"""
chr_dict={
    "NC_000001.11" : "chr1",
    "NC_000002.12" : "chr2",
    "NC_000003.12" : "chr3",
    "NC_000004.12" : "chr4",
    "NC_000005.10" : "chr5",
    "NC_000006.12" : "chr6",
    "NC_000007.14" : "chr7",
    "NC_000008.11" : "chr8",
    "NC_000009.12" : "chr9",
    "NC_000010.11" : "chr10",
    "NC_000011.10" : "chr11",
    "NC_000012.12" : "chr12",
    "NC_000013.11" : "chr13",
    "NC_000014.9" : "chr14",
    "NC_000015.10" : "chr15",
    "NC_000016.10" : "chr16",
    "NC_000017.11" : "chr17",
    "NC_000018.10" : "chr18",
    "NC_000019.10" : "chr19",
    "NC_000020.11" : "chr20",
    "NC_000021.9" : "chr21",
    "NC_000022.11" : "chr22",
    "NC_000023.11" : "chrX",
    "NC_000024.10" : "chrY"
    }


def load_deepsignal(filename):
    out = dict()
    pos_list = []
    #load data and for every line in input do...
    with open(filename) as filename_load:
        csv_reader = csv.reader(filename_load, delimiter='\t')
        for record in csv_reader:
            if len(record) != 0:

                #check if input is in valid chromosome format (no wrongly matched + mitochondrial DNA)
                if record[0] in chr_dict.keys():

                    #get information
                    chromosome = chr_dict[record[0]]
                    key=make_key(chromosome, record[1]) #record von 3 richtig? oder 1???
                    num_reads = record[8]
                    num_methylated = record[6]


                    #write information to the dictionary
                    #if new position: make new entry
                    if key not in out.keys():
                        pos_list.append(key)
                        out[key] = MethylationStats(chromosome, record[1], num_reads, num_methylated)

                    #if position already in there: sum up values
                    else:
                        out[key].accumulate(num_reads, num_methylated)

    return out, pos_list

"""
    Format of tombo_valid_coverage_file.wig

        variableStep chrom=NC_000001.11 span=1
        1672095 1
        ....
"""

def load_tombo(filename):
    out = dict()
    chromosome = ""
    suitable_chromosome = False
    pos_list = []

    #load data and for every line in input do...
    with open(filename) as filename_load:
        csv_reader = csv.reader(filename_load, delimiter=' ')
        for record in csv_reader:
            if len(record) != 0:

                #if new variable step -> get chromosome
                if record[0] == "variableStep":
                    split=record[1].split("=")
                    if split[1] in chr_dict.keys():
                        suitable_chromosome = True
                        chromosome = chr_dict[split[1]]
                    else:
                        suitable_chromosome = False

                #if position
                else:
                    #if the chromosome has been correctly identified
                    if suitable_chromosome == True:
                        position = to_int(record[0])-1
                        key=make_key(chromosome, position)
                        num_reads = 1
                        num_methylated = record[1]
                        #chromosome, position, num_reads, num_methylated
                        if key not in out.keys():
                            pos_list.append(key)
                            out[key] = MethylationStats(chromosome, position, num_reads, num_methylated)

                        #if position already in there: sum up values
                        else:
                            out[key].accumulate(num_reads, num_methylated)

    return out, pos_list


"""

    nanopolish format:

        chromosome | start | end | num_motifs_in_group | called_sites |
        called_sites_methylated | methylated_frequency | group_sequence

"""

def load_nanopolish(filename):
    out = dict()
    pos_list = []
    #load data and for every line in input do...
    with open(filename) as filename_load:
        csv_reader = csv.reader(filename_load, delimiter='\t')
        for record in csv_reader:
            if len(record) != 0:

                chromosome = record[0]
                position = record[1]-1
                key=make_key(chromosome, position)
                num_reads = record[4]
                num_methylated = record[5]

                if key not in out.keys():
                    pos_list.append(key)
                    out[key] = MethylationStats(chromosome, position, num_reads, num_methylated)

                #if position already in there: sum up values
                else:
                    out[key].accumulate(num_reads, num_methylated)

    return out, pos_list

"""
WGBS format:
Reference chromosome or scaffold | 0
Start position in chromosome | 1
End position in chromosome | 2
Name of item | 3
Score from 0-1000. Capped number of reads | 4
Strandedness, plus (+), minus (-), or unknown (.) | 5
Start of where display should be thick (start codon) | 6
End of where display should be thick (stop codon) | 7
Color value (RGB) | 8
Coverage, or number of reads | 9
Percentage of reads that show methylation at this position in the genome | 10
"""
chromosome_sammler=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
         "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
         "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

def load_bisulfite(filename):
    out = dict()
    pos_list = []
    #load data and for every line in input do...
    with open(filename) as filename_load:
        csv_reader = csv.reader(filename_load, delimiter='\t')
        for record in csv_reader:
            if len(record) != 0 and record[0] in chromosome_sammler:

                #get information from current row
                chromosome = record[0]
                position = record[1]
                key=make_key(chromosome, position)
                num_reads = float(record[9])
                methylation_frequency = float(record[10]) / 100
                num_methylated = proper_round(num_reads * methylation_frequency)
                if num_reads < num_methylated:
                    num_reads = num_methylated

                if key not in out.keys() and num_reads != 0:
                    pos_list.append(key)
                    out[key] = MethylationStats(chromosome, position, num_reads, num_methylated)

                #if position already in there: sum up values
                elif num_reads != 0:
                    out[key].accumulate(num_reads, num_methylated)

    return out, pos_list

"""
    needed methrix format:
        chromosome | position | methylated_reads | all_reads
"""
def write_methrix_output(file_dict, output_log_file):
    with open(output_log_file, 'w+', newline='') as log_file:
        tsv_writer = csv.writer(log_file, delimiter='\t')
        for pos, stats in file_dict.items():
            printer=[stats.chromosome(), stats.position(), stats.methylated_reads(), stats.unmethylated()]
            tsv_writer.writerow(printer)


#check_loading_style
def check_loading_style(filename):
    fstr = str(filename)
    if "deepsignal" in fstr:
        print("loading deepsignal")
        return load_deepsignal(filename)
    elif "tombo" in fstr:
        print("loading tombo")
        return load_tombo(filename)
    elif "nanopolish" in fstr:
        print("loading nanopolish")
        return load_nanopolish(filename)
    elif "bisulfite" in fstr:
        print("loading bisulfite")
        return load_bisulfite(filename)
    else:
        print("couldn`t identify filetype: need tombo, deepsignal, nanopolish or bisulfite files")


def intersection_of_2lists(lst1, lst2):
    keys1=set(lst1)
    keys2=set(lst2)
    intersect=keys1 & keys2
    return intersect

def write_frequency(set1, set2, pos_list, out_file):
    with open(out_file, 'w+', newline='') as log_file:
        tsv_writer = csv.writer(log_file, delimiter='\t')
        for pos in pos_list:
            d1=set1[pos]
            d2=set2[pos]
            d1_m = d1.methylation_frequency()
            d2_m = d2.methylation_frequency()
            printer=[pos, d1_m, d2_m]
            tsv_writer.writerow(printer)


def main():
    parser = argparse.ArgumentParser(description='read frequency file and write into frequency comparison format')
    parser.add_argument('--input_file', '-i', action="store", type=str, required=True,
                        help='the input file')
    parser.add_argument('--result_file', '-r', action="store", type=str, required=True,
                        help='the file path to save the result')
    args = parser.parse_args()


    bs_path = "/icgc/dkfzlsdf/analysis/C010/mayerma/methylation_analysis/output/all_bisulfite_frequency.tsv"
    file1, pos_1 = check_loading_style(bs_path)
    file2, pos_2 = check_loading_style(args.input_file)

    pos_list = intersection_of_2lists(pos_1, pos_2)
    write_frequency(file1, file2, pos_list, args.result_file)



    #write_methrix_output(file, args.result_file)

if __name__ == '__main__':
    sys.exit(main())
