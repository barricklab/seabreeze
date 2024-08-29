#!/usr/bin/env python3
# Author: Ira Zibbu
# Last update: 2024-04-16
# This script accepts a genomediff file, bin size and genome size to generate a 1-0 binary assignments for non-overlapping bins over the genome

''' imports '''

import os
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import argparse


''' fetch arguments from the command line '''

parser = argparse.ArgumentParser(description='parallel_rearrangement_stack.py')
parser.add_argument('--gd', help='path to GenomeDiff file')
parser.add_argument('--output', help='name of output csv file')
parser.add_argument('--bin', help='size of non-overlapping bins')
parser.add_argument('--size', help='total length of the genome')
parser.add_argument('--deletion',action='store_true', help='Only report deletions')
parser.add_argument('--duplication',action='store_true', help='Only report duplications')

def read_genome_diff(gd):

    ''' read the genome diff file into a dataframe '''
    data = []

    #  excluding lines starting with '#=' these are header lines
    with open(gd, 'r') as file:
        for line in file:
            if not line.startswith('#='):
                data.append(line.strip().split('\t'))

    df_gd = pd.DataFrame(data)
    return df_gd

def generate_bins(bin,size):

    ''' returns a nx2 dataframe of n bins of size 'bin' spread over genome of 'size' '''
    size = int(size)
    bin = int(bin)
    num_bins = size // bin + 1 # number of bins
    bins_list = list(range(0, size + 1, bin)) # list of bins
    
    df_stack = pd.DataFrame({'bin': bins_list, 'count': [0] * num_bins})
    return df_stack

def annotate_bins(df_gd, df_stack, deletion, duplication,bin):
    
    ''' annotate the bins as 1-0 depending on if it overlaps with a deletion '''
    bin = int(bin)
    if deletion: #only create stack for deletions
        for index, row in df_gd.iterrows():
            mutation_type = row[0]
            start_point = int(row[4])
            length = int(row[5])
        
            start_bin = start_point // bin
            end_bin = (start_point + length) // bin
            
            # Update df_stack accordingly
            if mutation_type == 'DEL':
                df_stack.loc[start_bin:end_bin, 'count'] = 1
                
    return df_stack

def main(gd,output,bin,size,deletion,duplication):
    df_gd = read_genome_diff(gd)
    df_stack = generate_bins(bin,size)
    df_stack = annotate_bins(df_gd,df_stack,deletion,duplication,bin)
    df_stack.to_csv(output, sep='\t', index=False)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.gd,args.output,args.bin,args.size,args.deletion,args.duplication)
