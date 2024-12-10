#!/usr/bin/env python3

"""
Unit test for plotsr
"""

__author__ = "Ira Zibbu"
__version__ = "0.0.1"

import os
import subprocess
import shutil
import bin.scripts.align_and_visual_timeseries as avt

def test_plotsr_genomes_tables():

    ''' test to see that the tsv file needed by plotsr is made '''

    subject_file="test/data/Anc-_0gen_REL606_truncated.fasta"
    query_file="test/data/Ara+1_10000gen_4530A_truncated.fasta"
    genomes_tsv_outfile="file1.genomes.tsv"

    fasta=[subject_file,query_file]
    names=["Anc-_0gen_REL606_truncated","Ara+1_10000gen_4530A_truncated"]
    count=1

    avt.make_genomes_file(fasta,names,count)

    assert os.path.isfile(genomes_tsv_outfile)
    assert os.path.getsize(genomes_tsv_outfile) > 0


def test_plotsr_plots():

    syri="test/data/file1syri.out"
    genomes="file1.genomes.tsv"
    count=1
    batch=True
    length=10
    plotsr="bin/scripts/plotsr/plotsr-bin"

    plot_pdf_file="file1_plotsr.pdf"
    genomes_tsv_outfile="file1.genomes.tsv"


    avt.makeplots(syri, genomes, count, batch,length,plotsr)

    assert os.path.isfile(plot_pdf_file)
    assert os.path.getsize(plot_pdf_file) > 0

    os.remove(plot_pdf_file)
    os.remove("plotsr.log")
    os.remove(genomes_tsv_outfile)
