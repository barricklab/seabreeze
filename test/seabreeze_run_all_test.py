"""
Master test to see if the seabreeze 'run_all' workflow pipelines run
"""

import subprocess
import pytest
import os
import shutil

def test_run_run_all_masked():

    command = ["seabreeze", "run", "run_all", "--dir","test/seabreeze_run_test","--masked","--assembly", "test/seabreeze_run_test/test_genomes/Ara+1_10000gen_4530A_truncated.fasta", "--ancestor", "test/seabreeze_run_test/test_genomes/Anc-_0gen_REL606_truncated.fasta", "--ori","AGGATGCTTTACCCAATATCAGCGAT","--dif","CTGGTCGGCAGAATGAGCAATCGCCA"]
    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    files_to_remove = ["test/seabreeze_run_test/03_reindex_genomes","test/seabreeze_run_test/04_rename_genome","test/seabreeze_run_test/05_isescan_tables","test/seabreeze_run_test/06_nucmer_alignment","test/seabreeze_run_test/07_syri_output","test/seabreeze_run_test/08_reindex_genome_oric","test/seabreeze_run_test/09_annotated_genomes","test/seabreeze_run_test/11_annotated_boundaries","test/seabreeze_run_test/12_genome_diff_tables"]

    for file in files_to_remove:
        shutil.rmtree(file)

def test_run_run_all_unmasked():

    command = ["seabreeze", "run", "run_all","--dir","test/seabreeze_run_test", "--masked","--assembly", "test/seabreeze_run_test/test_genomes/Ara+1_10000gen_4530A_truncated.fasta", "--ancestor", "test/seabreeze_run_test/test_genomes/Anc-_0gen_REL606_truncated.fasta", "--ori","AGGATGCTTTACCCAATATCAGCGAT","--dif","CTGGTCGGCAGAATGAGCAATCGCCA"]
    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    files_to_remove = ["test/seabreeze_run_test/03_reindex_genomes","test/seabreeze_run_test/04_rename_genome","test/seabreeze_run_test/05_isescan_tables","test/seabreeze_run_test/06_nucmer_alignment","test/seabreeze_run_test/07_syri_output","test/seabreeze_run_test/08_reindex_genome_oric","test/seabreeze_run_test/09_annotated_genomes","test/seabreeze_run_test/11_annotated_boundaries","test/seabreeze_run_test/12_genome_diff_tables"]

    for file in files_to_remove:
        shutil.rmtree(file)
