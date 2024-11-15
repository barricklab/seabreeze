#!/usr/bin/env python3
'''
Unit tests for reindex_assembly.py
'''

__author__ = "Ira Zibbu"
__version__ = "0.0.1"

'''
imports
'''

from bin.scripts.reindex_assembly import reindex_fasta
from Bio.Seq import Seq
from Bio import SeqIO

ltee_test_data="test/data/Ara-5_75000gen_A_ori.1.fasta"
ltee_test_data_output="test/data/Ara-5_75000gen_A.fasta"
reindex_bases="AGCTTTTCATTC"

def test_reindexing_ltee_genome():

    reindex_fasta(ltee_test_data,ltee_test_data_output,reindex_bases)
    fasta = SeqIO.read(ltee_test_data_output,"fasta")
    first_12_bases=fasta.seq[:12]
    assert first_12_bases=="AGCTTTTCATTC"
