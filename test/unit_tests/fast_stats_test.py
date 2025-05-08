#!/usr/bin/env python3

"""
Unit test for rename_contigs.py
"""

__author__ = "Ira Zibbu"
__version__ = "0.0.1"

from src.seabreeze.fasta_stats import get_genome_sizes
import os
from Bio import SeqIO

def test_rename_contigs():
    yersinia_KIM="test/data/yersinia_pestis_KIM.fasta"
    yersinia_nepal514="test/data/yersinia_pestis_nepal514.fasta"
    anc0gen="test/data/Anc-_0gen_REL606.fasta"
    araplus1="test/data/Ara+1_10000gen_4530A.fasta"

    yersinia_KIM_fasta=SeqIO.read(yersinia_KIM,"fasta")
    yersinia_nepal514_fasta=SeqIO.read(yersinia_nepal514,"fasta")
    anc0gen_fasta=SeqIO.read(anc0gen,"fasta")
    araplus1_fasta=SeqIO.read(araplus1,"fasta")

    genome_stats_yersinia=get_genome_sizes(yersinia_nepal514_fasta,yersinia_KIM_fasta,)
    genome_stats_ltee=get_genome_sizes(araplus1_fasta,anc0gen_fasta)

    assert genome_stats_yersinia['size_assembly']==4534590
    assert genome_stats_yersinia['size_ancestor']==4600755
    assert genome_stats_yersinia['difference']==-66165

    assert genome_stats_ltee['size_assembly']==4650661
    assert genome_stats_ltee['size_ancestor']==4629812
    assert genome_stats_ltee['difference']==20849
