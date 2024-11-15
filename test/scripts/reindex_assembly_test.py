#!/usr/bin/env python3
"""
Test for reindex_assembly.py
"""


__author__ = "Ira Zibbu"
__version__ = "0.0.1"

''' imports '''
import argparse
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description='run test on output of reindex_assembly.py')
parser.add_argument('--file', help='path to the FASTA file with the assembly reindexed')
parser.add_argument('--output', help='output file to write results of the test to')

def test(file):

    assembly = SeqIO.read(file, "fasta")
    first_12_bases=assembly.seq[:12]

    # These conditions are true from prior runs and is the expected result
    if "yersinia" in file:
        assert first_12_bases=="TCGCGCGATCTT",f"{file} does not begin at the correct bases. Reindexing failed"
        return True


    if "gen" in file:
        assert first_12_bases=="AGCTTTTCATTC",f"{file} does not begin at the correct bases. Reindexing failed"
        return True

def main(file,output):

    flag = test(file)
    if flag:
        with open(output, "w") as file:
            file.write("test successful")


if __name__ == '__main__':
    args = parser.parse_args()
    main(args.file,args.output)
