#!/usr/bin/env python3

"""
Unit test for replichore_arms_analyse.py
"""

__author__ = "Ira Zibbu"
__version__ = "0.0.1"

from src.seabreeze.replichore_arms_analyse import coords

def test_finding_coords_REL606():

    file = "test/data/Anc-_0gen_REL606.fasta"
    ori = "GGATCCTGGGTATTAAAA"
    dif= "TCTTCCTTGGTTTATATT"

    REL606_coords = coords(file,ori,dif)

    assert REL606_coords["oric_start"]==3886082
    assert REL606_coords["dif_start"]==1567362

def test_finding_coords_REL606_evolved_1():

    file = "test/data/REL606_evolved_1.fasta"
    ori = "GGATCCTGGGTATTAAAA"
    dif= "TCTTCCTTGGTTTATATT"

    REL606_evoled_1_coords = coords(file,ori,dif)

    assert REL606_evoled_1_coords["oric_start"]==3876626
    assert REL606_evoled_1_coords["dif_start"]==1232116
