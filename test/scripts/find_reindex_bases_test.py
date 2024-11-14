#!/usr/bin/env python3
"""
Test for find_reindex_bases.py
"""

__author__ = "Ira Zibbu"
__version__ = "0.0.3"

''' imports '''
import argparse

parser = argparse.ArgumentParser(description='run test on output of find_reindex_bases rule')
parser.add_argument('--file', help='txt file with bases output by find_reindex_bases.py')
parser.add_argument('--output', help='output file to write to')

def test(file):

    with open(file, "r") as file_object:
        for line in file_object:
            bases=line.strip()

        # These conditions are true from prior runs and is the expected result
        if "yersinia" in file:
            assert bases=="TCGCGCGATCTT",f"Reindex bases were not found in {file}"
            return True

        if "gen" in file:
            assert bases=="AGCTTTTCATTC",f"Reindex bases were not found in {file}"
            return True

def main(file,output):

    flag = test(file)
    if flag:
        with open(output, "w") as file:
            file.write("test successful")

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.file,args.output)
