#!/usr/bin/env python3

"""
Unit test for syri
"""

__author__ = "Ira Zibbu"
__version__ = "0.0.1"

import os
import subprocess
import shutil
import bin.scripts.align_and_visual_timeseries as avt

def test_variant_calling():

    ''' test to see that the final delta coords file is generated and is not empty'''


    delta_filtered="test/data/file1.filtered.delta"
    delta_coords="test/data/file1.filtered.coords"
    subject_file="test/data/Anc-_0gen_REL606_truncated.fasta"
    query_file="test/data/Ara+1_10000gen_4530A_truncated.fasta"

    syri_outfile="file1syri.out"

    avt.call_variants(delta_coords, delta_filtered, subject_file, query_file, 1)

    assert os.path.isfile(syri_outfile)
    assert os.path.getsize(syri_outfile) > 0
    os.remove(syri_outfile)
