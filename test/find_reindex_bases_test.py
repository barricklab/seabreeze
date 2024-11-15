#!/usr/bin/env python3
'''
Unit tests for find_reindex_bases.py
'''

__author__ = "Ira Zibbu"
__version__ = "0.0.1"

'''
imports
'''

from bin.scripts.find_reindex_bases import find_smallest_unique_k_mer
from Bio.Seq import Seq
from Bio import SeqIO

yersinia_test_data="test/data/yersinia_pestis_KIM.fasta"
ltee_test_data="test/data/Anc-_0gen_REL606.fasta"

def test_check_unique_kmer_for_yersinia():

        fasta = SeqIO.read(yersinia_test_data, "fasta")
        k_mer = find_smallest_unique_k_mer(fasta)
        assert k_mer=="TCGCGCGATCTT"

def test_check_unique_kmer_for_ltee():

        fasta = SeqIO.read(ltee_test_data, "fasta")
        k_mer = find_smallest_unique_k_mer(fasta)
        assert k_mer=="AGCTTTTCATTC"
