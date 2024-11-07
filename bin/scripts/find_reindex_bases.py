#!/usr/bin/env python3
# This script accepts a query and a subject sequence, and finds the first k-mer from the subject to reindex the query sequence to. The query sequence is rotated by a separate script

__author__ = "Ira Zibbu"
__version__ = "0.0.2"

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
parser.add_argument('--output', help='output txt files with bases to reindex to')

def load_test_fasta_files(subject_path, query_path):
    ''' Accepts the paths to the subject and query, checks if they are in the correct format and loads them as a Seq objects '''

    if not(os.path.exists(subject_path)):
        raise IOError(f"{subject_path} subject file does not exist")

    if not(os.path.exists(query_path)):
        raise IOError(f"{query_path} subject file does not exist")

    subject = None
    query = None

    try:
        subject = SeqIO.read(subject_path, "fasta")
    except Exception as e:
        print(f"Error parsing subject file: {subject_path}. Exception: {e}")
        raise

    try:
        query = SeqIO.read(query_path, "fasta")
    except Exception as e:
        print(f"Error parsing query file: {query_path}. Exception: {e}")
        raise

    if subject is None or query is None:
        raise ValueError("Both subject and query files must contain valid FASTA records.")

    return subject, query

    # try:
    #     subject = SeqIO.read(subject_path, "fasta")
    #     query = SeqIO.read(query_path, "fasta")
    # except:
    #     print("Error parsing input files. Check that files are in FASTA format and have only one non-empty record")

    # return subject, query

def check_k_mer_unique(fasta,k_mer):

    ''' Accepts a SeqRecord and a k-mer, and checks if that k-mer is unique in the sequence '''

    frequency_of_k_mer = fasta.count(k_mer)
    if frequency_of_k_mer == 1:
        #print(f" {k_mer} is unique")
        return True

def find_smallest_unique_k_mer(fasta):

    ''' Finds the smallest k-mer starting at from the 1st base of the sequence that is unique '''

    min_k_value = 10 # this value is arbitrary but comes from my experience with E. coli genomes
    max_k_value = 100 # see above comment

    k_mer_found = False
    for n in range(min_k_value,max_k_value):
        k_mer = fasta[0:n].seq
        if check_k_mer_unique(fasta,k_mer):
            k_mer_found = True
            return k_mer

    if not(k_mer_found):
        raise Exception("No unique k-mer found. Please manually specify a sequence")

def main(subject_path, query_path,output):
    subject, query = load_test_fasta_files(subject_path,query_path)
    k_mer = find_smallest_unique_k_mer(subject)
    if check_k_mer_unique(query,k_mer):
        print(f"The smallest k-mer unique in both subject and query is {len(k_mer)} with sequence {k_mer}")
    with open(output, "w") as file:
        file.write(str(k_mer))

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.subject, args.query,args.output)
