#!/usr/bin/env python3

"""
Unit test for rename_contigs.py
"""

__author__ = "Ira Zibbu"
__version__ = "0.0.1"

from bin.scripts.rename_contigs import rename_contigs
from Bio import SeqIO
import os

def test_rename_contigs():
    file="test/data/Anc-_0gen_REL606.fasta"
    output_file="genome.temp"
    name="genome"

    rename_contigs(file,name,output_file)
    fasta = SeqIO.read(output_file, "fasta")
    assert fasta.id=="genome"

    os.remove(output_file)
