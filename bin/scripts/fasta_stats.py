#!/usr/bin/env python3
# Author: Ira Zibbu
# this script outputs the number of contigs, the length of each contig, the length of the ancestor and the size difference

''' imports '''

import os
from collections import deque
import pandas as pd
import numpy as np
import argparse

''' fetch arguments '''

parser = argparse.ArgumentParser(description='fasta_stats.py, a script to calculate the number of contigs and length of contigs of fasta files, and size difference relative to ancestor')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')
parser.add_argument('--file', '-f', action='store_true', help='Enable single file mode. Default is folder mode')
parser.add_argument('--fasta', help='FASTA file. Use with the --file option')
parser.add_argument('--folder', help='Folder of FASTA files. Do not use with --file option')
parser.add_argument('--ancestor', help='FASTA file of the ancestor, relative to which genome size change is calculated')
parser.add_argument('--output', help='Output filename for table with contig sizes')
parser.add_argument('--stats', help='Output filename for table with genome size difference. Use with the ancestor option')


def get_fasta_names(folder_path):

    ''' fetch names of fasta files in the folder '''

    fasta_names = []
    for entry in os.scandir(folder_path):
        if entry.is_file() and entry.name.endswith(".fasta") and len(entry.name) > len(".fasta"):
            fasta_names.append(entry.name)
    return fasta_names

def fasta_length(name):

    ''' return the names and lengths of contigs from name '''

    fasta = [] # a list of deques, where each deque corresponds to one contig, and the elements of each deque contains the lines of the contig
    contig_lengths = [] # a list of the lengths of the contigs
    contig_count=0 # number of contigs 
    contig_names=[] # a list of the names of the contigs

    with open(name, "r") as file:
        for line in file:
            line = line.strip()
            if line[0] == ">":
                contig_count+=1
                contig_names.append(line)
                fasta.append(deque()) # make a deque for this contig
                contig_lengths.append(0)
            else:
                data = line
                fasta[(contig_count-1)].append(data) # append the current line of the fasta file into the deque of this contig
        
        for count, element in enumerate(fasta):
            seq=''
            seq=''.join(element) # combine the elements of the deque into a single sequence i.e seq is a string of the sequence of the contig
            contig_lengths[count]=len(seq)

    return contig_names, contig_lengths

def contig_count(fasta_names):

    ''' determines the maximum number of contigs among all fasta files provivded '''

    max_contig_count=0
    for name in fasta_names:
        with open(name, "r") as file:
            contig_count=0
            for line in file:
                line = line.strip()
                if line[0] == ">":
                    contig_count+=1
            max_contig_count=max(max_contig_count,contig_count)  
    return max_contig_count    
        
def analyse_genomes(df, ancestor_filename, stats):

    ''' accepts a dataframe, name of ancestor to calculate genome size differences '''

    ancestor=df[df['File'].str.contains(ancestor_filename)]
    ancestor_dict = ancestor.iloc[0].to_dict()
    ancestor_length=ancestor_dict['contig_1']
    col_names=['clone', 'length', 'difference', 'change']
    df_genomes=pd.DataFrame(np.nan, index=range(len(df)), columns=col_names)
    df_genomes['clone']=df['File']
    df_genomes['clone'] = df_genomes['clone'].str.replace('.fasta', '')
    df_genomes['length']=df['contig_1']
    for idx in range(len(df_genomes)):
        difference=int(df_genomes.loc[idx,'length'])-int(ancestor_length)
        df_genomes.loc[idx,'difference']=difference
        change=(difference/int(ancestor_length))*100
        df_genomes.loc[idx,'change']=change
    print(df_genomes.head)

    stats_dir = os.path.dirname(stats)

    try:
        # Save the DataFrame to CSV
        df_genomes.to_csv(stats, sep='\t', index=False)
    except OSError as e:
        raise OSError(f"Cannot save file into a non-existent directory: '{stats_dir}'. Error: {str(e)}")

    #df_genomes.to_csv(stats, sep='\t', index=False)


def main(file, fasta, folder, ancestor, output, stats):

    name='' 
    print("starting task")

    if not file: # folder mode

        fasta_names = get_fasta_names(folder)
        wd = os.getcwd()
        os.chdir(folder)
        max_contig_count=contig_count(fasta_names)
        col_names=['File'] 

        for idx in range(max_contig_count): # generate column names as contig_1, comntig_2.. contig_n where n is the mac contig count
            idx+=1
            name="contig_"+str(idx)
            col_names.append(name)
            idx-=1

        df_contigs=pd.DataFrame(np.nan, index=range(len(fasta_names)), columns=col_names)
        
        row_idx=0
        for name in fasta_names:
            contig_names, contig_lengths=fasta_length(name)
            df_contigs.loc[row_idx,'File']=name
            col_idx=1
            for num in contig_lengths: # add the lengths of the contigs in the row for name
                contig_col=col_names[col_idx]
                df_contigs.loc[row_idx,contig_col]=str(num)
                col_idx+=1
            row_idx+=1
        
        os.chdir(wd) # go back to the directory 
        
        df_contigs.to_csv(output, sep='\t', index=False)

        output_dir = os.path.dirname(output)

        try:
            # Save the DataFrame to CSV
            df_contigs.to_csv(output, sep='\t', index=False)
        except OSError as e:
            raise OSError(f"Cannot save file into a non-existent directory: '{output_dir}'. Error: {str(e)}")

        analyse_genomes(df_contigs, ancestor, stats)

    if file: # file mode
    
        contig_names, contig_lengths=fasta_length(fasta)
        for count, element in enumerate(contig_names): # loop to print out contig names and contig lengths
            print(element, " : ", contig_lengths[count])
    
        print("Task finished")


if __name__ == '__main__':
    args = parser.parse_args()
    main(args.file, args.fasta, args.folder, args.ancestor, args.output, args.stats)