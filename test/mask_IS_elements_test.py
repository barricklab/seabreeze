#!/usr/bin/env python3
"""
unit test for mask_IS_elements_test.py
"""
__author__ = "Ira Zibbu"
__version__ = "0.0.1"

''' imports '''

from Bio.Seq import Seq
from Bio import SeqIO
import bin.scripts.mask_IS_elements as mask_IS_elements

def test_masking_bases():
    fasta_file="test/data/Anc-_0gen_REL606.fasta"
    #is_table_file="test/data/Anc-0gen_REL606.csv"
    fasta=SeqIO.read(fasta_file,"fasta")

    # picked two random IS elements in this genome
    start1=15386
    stop1=16728

    start2=149489
    stop2=150576

    masked_fasta_1=mask_IS_elements.mask_bases(fasta,start1,stop1)
    predicted_bases='N' * (stop1 - start1)
    tmp=str(masked_fasta_1.seq[start1:stop1]) # just a string representing the sequence of masked bases

    assert tmp==predicted_bases

    masked_fasta_1=mask_IS_elements.mask_bases(fasta,start2,stop2)
    predicted_bases='N' * (stop1 - start1)
    tmp=str(masked_fasta_1.seq[start1:stop1]) # just a string representing the sequence of masked bases

    assert tmp==predicted_bases
