#!/usr/bin/env python3
"""
Take the output table from ISEscan, and mask those bases in the supplied genome
"""
__author__ = "Ira Zibbu"
__version__ = "0.0.1"

''' imports '''
import os
import argparse
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO

''' fetch arguments '''
parser = argparse.ArgumentParser(description='mask_IS_elements.py, a script to mask IS elements with Ns')
parser.add_argument('--genome', help='path to the FASTA file')
parser.add_argument('--is_table', help='path to the ISEscan output csv file')
parser.add_argument('--output', help="output file to write modified genome to")

def read_is_table(is_table_file):

    """ Accept path to csv file and return a dataframe """

    is_table=pd.read_csv(is_table_file)

    if is_table.shape[1] != 24:
        raise IOError(f"{is_table_file} has an incorrect number of columns. Please check format")

    return is_table

def read_fasta(genome):

    """ Read fasta file and return seq record object """

    if not(os.path.exists(genome)):
        raise IOError(f"{genome} file does not exist")

    try:
        fasta = SeqIO.read(genome, "fasta")
    except Exception as e:
        print(f"Error parsing FASTA file: {genome}. Exception: {e}")
        raise

    return fasta

def mask_bases(fasta,start,stop):

    """ Accept a seqrecord and start/stop site and replace those bases with NNNs"""

    modified_sequence = fasta.seq[:start] + 'N' * (stop - start + 1) + fasta.seq[stop:]
    fasta.seq=modified_sequence
    return fasta

def main(genome,is_table_file,output):

    is_table=read_is_table(is_table_file)
    fasta=read_fasta(genome)

    # interate over the lines of the ISEscan output and mask each IS element

    modified_fasta=fasta
    for index, row in is_table.iterrows():
        start=row["isBegin"]
        stop=row["isEnd"]
        is_name=row["cluster"]
        modified_fasta=mask_bases(modified_fasta,start,stop)
        print(f"Masked IS {is_name} from {start} to {stop}")

    SeqIO.write(modified_fasta, output, "fasta")

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.genome, args.is_table,args.output)
