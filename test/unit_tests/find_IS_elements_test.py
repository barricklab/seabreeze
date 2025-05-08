#!/usr/bin/env python3

"""
Unit test for isescan
"""

__author__ = "Ira Zibbu"
__version__ = "0.0.1"

import os
import subprocess
import shutil

def test_isecan():
    file="test/data/Anc-_0gen_REL606_truncated.fasta"
    output_csv="test/data/Anc-_0gen_REL606_truncated/data/Anc-_0gen_REL606_truncated.fasta.csv"
    output_folder = "test/data/Anc-_0gen_REL606_truncated"

    if os.path.isfile(output_csv):
        shutil.rmtree(output_folder)

    command=["isescan.py", "--seqfile",file, "--output",output_folder]
    subprocess.run(command, check=True)
    assert os.path.isfile(output_csv)
    shutil.rmtree(output_folder)
