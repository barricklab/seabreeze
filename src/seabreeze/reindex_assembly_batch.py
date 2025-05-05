#!/usr/bin/env python3
"""
Accept a folder of fasta files, a data.csv and ori_dif_sequences.csv file
And reindex each assembly to its own ancestor
"""

__author__ = "Ira Zibbu"
__version__ = "0.1.0"

""" imports """
import os
import pandas as pd
import numpy as np
import argparse
from replichore_arms_analyse import get_sequence
from reindex_assembly import reindex_fasta
from fasta_stats import read_data


""" arguments """
parser = argparse.ArgumentParser(description='what does this script to?')
parser.add_argument('--folder', help='folder of fasta files')
parser.add_argument('--data', help='path to the data.csv file')
parser.add_argument('--sequences', help='path to the csv file with the sequences for the ori and dif loci')
parser.add_argument('--output', help='Output folder to put reindexed genomes in')


def main(folder,data,sequences,output):

    data_df=read_data(data) # table of genome assemblies and their ancestor
    sequences_df=read_data(sequences) # table of ancestors and their ori and dif sequences
    assembly_to_ancestor_dict = dict(zip(data_df['assembly'], data_df['ancestor']))
    assemblies=data_df['assembly'].tolist() # a list of all of the assemblies specified in the data.csv file

    for name in assemblies:

        file_name=f"{name}.fasta" # input file name
        file_in=os.path.join(folder,file_name)
        file_name=f"{name}.fasta" # output file name
        file_out=os.path.join(output,file_name)
        oric=get_sequence(sequences_df, assembly_to_ancestor_dict[name],'ori')
        reindex_fasta(file_in, file_out, oric)

if __name__ == "__main__":
	args = parser.parse_args()
	main(args.folder,args.data,args.sequences,args.output)
