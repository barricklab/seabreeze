#!/usr/bin/env python3
# Author: Ira Zibbu
# this script outputs the number of contigs, the length of each contig, the length of the ancestor and the size difference

__author__ = "Ira Zibbu"
__version__ = "0.1.0"

''' imports '''

import os
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
import argparse

''' fetch arguments '''

parser = argparse.ArgumentParser(description='fasta_stats.py, a script to calculate the number of contigs and length of contigs of fasta files, and size difference relative to ancestor')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')
parser.add_argument('--folder', help='Folder of FASTA files.')
parser.add_argument('--output', help='Output filename for table with genome size difference')
parser.add_argument('--data',help='data.csv file with subject and queries for pairwise comparisons ')

def read_data(data):

    ''' Accept the path to the data.csv file an return as a pandas dataframe '''

    try:
        path_to_csv=data
        data = pd.read_csv(path_to_csv)
    except Exception as e:
        print(f"Error parsing file: {data}. Exception: {e}")
        raise

    return data

def get_genome_sizes(assembly,ancestor):

    ''' Accept two seq record objects, and return the sizes, size difference and percent change in a genome_statsionary '''

    genome_stats = {'size_assembly':0,'size_ancestor':0,'difference':0,'percent_change':0}
    genome_stats['size_assembly'] = len(assembly.seq)
    genome_stats['size_ancestor'] = len(ancestor.seq)
    genome_stats['difference'] = genome_stats['size_ancestor']-genome_stats['size_assembly']
    genome_stats['percent_change'] = genome_stats['difference']/genome_stats['size_ancestor']

    return genome_stats

def main(verbose,folder,output,data):

    from find_reindex_bases import load_test_fasta_files
    # this import statement was moved to avoid missing imports during testing

    data = read_data(data)
    columns_for_table = ['assembly','ancestor','size_assembly','size_ancestor','difference','percent_change']
    genome_size_table = pd.DataFrame(columns=columns_for_table)

    # loop over data
    for index, row in data.iterrows():
        assembly_file=f"{row['assembly']}.fasta"
        ancestor_file=f"{row['ancestor']}.fasta"

        assembly_file_path=os.path.join(folder,assembly_file)
        ancestor_file_path=os.path.join(folder,ancestor_file)

        assembly, ancestor = load_test_fasta_files(assembly_file_path,ancestor_file_path)
        genome_stats = get_genome_sizes(assembly,ancestor)
        genome_size_table.loc[index,'assembly']=row['assembly']
        genome_size_table.loc[index,'ancestor']=row['ancestor']

        for key,value in genome_stats.items():
            genome_size_table.loc[index,key]=value

    try:
        genome_size_table.to_csv(output,index=False)
    except Exception as e:
        print(f"Error writing file {output}. Exception: {e}")
        raise

    print(f"Genome size statistics computed. Results stored in {output}")

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.verbose, args.folder, args.output,args.data)
