#!/usr/bin/env python3
# Author: Ira Zibbu
# this script accepts the name of a folder with FASTA files or a FASTA file to rename the all contigs to a specified string


__author__ = "Ira Zibbu"
__version__ = "0.0.2"

''' imports '''

import os
import glob
from Bio import SeqIO
import argparse

''' fetch arguments '''

parser = argparse.ArgumentParser(description='rename_contigs.py, a script to rename all the contigs of FASTA files to specified string')
parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose mode')
# parser.add_argument('--file', '-f', action='store_true', help='Enable single file mode. Default is folder mode')
parser.add_argument('--file', help='FASTA file. Use with the --file option')
# parser.add_argument('--folder', help='Folder of FASTA files. Do not use with --file option')
parser.add_argument('--output', help='Name of output fasta file. Use with --file and --fasta options ONLY')
parser.add_argument('--name', help='String to rename contigs to')

# def get_fasta_names(folder):

#     ''' return names of all fasta files in the specified directory '''

#     fasta_names = []
#     for entry in os.scandir(folder):
#         if entry.is_file() and entry.name.endswith(".fasta") and len(entry.name) > len(".fasta"):
#             #print(entry.name)
#             fasta_names.append(entry.name)
#     return fasta_names

def load_test_fasta_files(file):
    ''' Accepts the path to the FASTA file, checks if they are in the correct format and loads them as Seq objects '''

    if not(os.path.exists(file)):
        raise IOError(f"{file} file does not exist")

    fasta = None

    try:
        fasta = SeqIO.read(file, "fasta")
    except Exception as e:
        print(f"Error parsing file: {file}. Exception: {e}")
        raise

    if fasta is None:
        raise ValueError(f"{file} must contain valid FASTA records.")

def rename_contigs(file,name,output):

    ''' rename all contigs in fasta_file to name '''

    modified_records = []

    # Read the FASTA file and modify the header lines
    with open(file, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            record.id = record.id if record.id is not None else "default_id"
            record.id = name
            record.description = ""
            modified_records.append(record)

    with open(output, "w") as fasta_file:
        SeqIO.write(modified_records, fasta_file, "fasta")
        print("Renamed contigs to ", name, " in ", output)

def main(file,name,output):
    load_test_fasta_files(file)
    rename_contigs(file, name, output)

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.file, args.name,args.output)
