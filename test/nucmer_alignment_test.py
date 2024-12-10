#!/usr/bin/env python3

"""
Unit test for nucmer
"""

__author__ = "Ira Zibbu"
__version__ = "0.0.1"

import os
import subprocess
import shutil
import bin.scripts.align_and_visual_timeseries as avt

def test_aligner():

    ''' test to see that the final delta coords file is generated and is not empty'''

    subject_file="test/data/Anc-_0gen_REL606_truncated.fasta"
    query_file="test/data/Ara+1_10000gen_4530A_truncated.fasta"

    delta_outfile="file1.delta"
    delta_filtered_outfile="file1.filtered.delta"
    delta_coords_outfile="file1.filtered.coords"


    if os.path.isfile(delta_coords_outfile):
        os.remove(delta_coords_outfile)


    if os.path.isfile(delta_outfile):
            os.remove(delta_outfile)


    if os.path.isfile(delta_filtered_outfile):
            os.remove(delta_filtered_outfile)

    avt.aligner(subject_file,query_file,1)

    assert os.path.isfile(delta_coords_outfile)
    assert os.path.getsize(delta_coords_outfile) > 0
    os.remove(delta_coords_outfile)
    os.remove(delta_outfile)
    os.remove(delta_filtered_outfile)
