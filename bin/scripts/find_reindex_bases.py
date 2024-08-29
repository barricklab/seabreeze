#!/usr/bin/env python3
# This script accepts a query and a subject sequence, and finds the first k-mer from the subject to reindex the query sequence to. The query sequence is rotated by a separate script

__author__ = "Ira Zibbu"
__version__ = "0.0.1"

''' imports '''
import os
import argparse
from re import sub
from Bio.Seq import Seq
from Bio import SeqIO

''' fetch arguments '''
parser = argparse.ArgumentParser(description='find_reindex_bases.py, a script to find the first k-mer to rotate the query sequence to be in the same ')
parser.add_argument('--subject', help='path to the subject FASTA file')
parser.add_argument('--query', help='path to the query FASTA file')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')


def load_test_fasta_files(subject_path, query_path):
    ''' Accepts the paths to the subject and query, checks if they are in the correct format and loads them as a Seq objects '''

    if not(os.path.exists(subject_path)):
        raise IOError(f"{subject_path} subject file does not exist")

    if not(os.path.exists(query_path)):
        raise IOError(f"{query_path} subject file does not exist")

    try:
        subject = SeqIO.read(subject_path, "fasta")
        query = SeqIO.read(query_path, "fasta")
    except:
        print("Error parsing input files. Check that files are in FASTA format and have only one non-empty record")

    return subject, query

def check_k_mer_unique(fasta,k_mer):
    frequency_of_k_mer = fasta.count(k_mer)
    if frequency_of_k_mer == 1:
        #print(f" {k_mer} is unique")
        return True

def find_smallest_unique_k_mer(fasta):

    ''' Finds the smallest k-mer starting at from the 1st base of the sequence that is unique '''

    min_k_value = 10
    max_k_value = 31

    k_mer_found = False
    for n in range(min_k_value,max_k_value):
        k_mer = fasta[0:n].seq
        if check_k_mer_unique(fasta,k_mer):
            k_mer_found = True
            return k_mer

    if not(k_mer_found):
        raise Exception("No unique k-mer found. Please manually specify a sequence")

def main(subject_path, query_path):
    subject, query = load_test_fasta_files(subject_path,query_path)
    k_mer = find_smallest_unique_k_mer(subject)
    if check_k_mer_unique(query,k_mer):
        print(f"The smallest k-mer unique in both subject and query is {len(k_mer)} with sequence {k_mer}")


if __name__ == '__main__':
    args = parser.parse_args()
    main(args.subject, args.query)
