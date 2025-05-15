"""
Master test to see if the seabreeze 'run_all' workflow pipelines run
"""

import subprocess
import pytest
import os
import shutil

def test_batch_run_all_masked():

    command = ["seabreeze", "run", "run_all","--masked","--data", "test/seabreeze_batch_test/data.csv", "--dir", "test/seabreeze_batch_test/", "--oridif", "test/seabreeze_batch_test/ori_dif_sequences.csv"]
    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    files_to_remove = ["test/seabreeze_batch_test/03_reindex_genomes","test/seabreeze_batch_test/04_rename_genome","test/seabreeze_batch_test/05_isescan_tables","test/seabreeze_batch_test/06_nucmer_alignment","test/seabreeze_batch_test/07_syri_output","test/seabreeze_batch_test/08_reindex_genome_oric","test/seabreeze_batch_test/09_annotated_genomes","test/seabreeze_batch_test/11_annotated_boundaries","test/seabreeze_batch_test/12_genome_diff_tables"]

    for file in files_to_remove:
        shutil.rmtree(file)

def test_batch_run_all_unmasked():

    command = ["seabreeze", "batch", "run_all", "--data", "test/seabreeze_batch_test/data.csv", "--dir", "test/seabreeze_batch_test/", "--oridif", "test/seabreeze_batch_test/ori_dif_sequences.csv"]
    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    files_to_remove = ["test/seabreeze_batch_test/03_reindex_genomes","test/seabreeze_batch_test/04_rename_genome","test/seabreeze_batch_test/05_isescan_tables","test/seabreeze_batch_test/06_nucmer_alignment","test/seabreeze_batch_test/07_syri_output","test/seabreeze_batch_test/08_reindex_genome_oric","test/seabreeze_batch_test/09_annotated_genomes","test/seabreeze_batch_test/11_annotated_boundaries","test/seabreeze_batch_test/12_genome_diff_tables"]

    for file in files_to_remove:
        shutil.rmtree(file)
