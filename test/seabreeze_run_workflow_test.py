"""
Master test to see if the seabreeze run workflow pipelines runs
"""

import subprocess
import pytest
import csv
import os
import shutil

def read_csv_as_dicts(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))

def test_run_analyse_genome_sizes():

    command = ["seabreeze", "run", "analyse_genome_sizes", "--dir","test/seabreeze_run_test","--assembly", "test/seabreeze_run_test/test_genomes/Ara+1_10000gen_4530A_truncated.fasta", "--ancestor", "test/seabreeze_run_test/test_genomes/Anc-_0gen_REL606_truncated.fasta", "--ori","AGGATGCTTTACCCAATATCAGCGAT","--dif","CTGGTCGGCAGAATGAGCAATCGCCA", "--conda-frontend", "conda"]
    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    # path_to_template_output="test/seabreeze_run_test/expected/genome_size_stats.csv"
    # path_to_generated_output="test/seabreeze_run_test/04_rename_genome/genome_size_stats.csv"

    # actual = read_csv_as_dicts(path_to_generated_output)
    # expected = read_csv_as_dicts(path_to_template_output)

    # assert actual == expected

def test_run_predict_IS_elements():

    command = ["seabreeze", "run", "predict_IS_elements","--dir","test/seabreeze_run_test", "--assembly", "test/seabreeze_run_test/test_genomes/Ara+1_10000gen_4530A_truncated.fasta", "--ancestor", "test/seabreeze_run_test/test_genomes/Anc-_0gen_REL606_truncated.fasta", "--ori","AGGATGCTTTACCCAATATCAGCGAT","--dif","CTGGTCGGCAGAATGAGCAATCGCCA", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    # path_to_template_output="test/seabreeze_run_test/expected/Anc-_0gen_REL606_truncated.csv"
    # path_to_generated_output="test/seabreeze_run_test/05_isescan_tables/Anc-_0gen_REL606_truncated.csv"

    # actual = read_csv_as_dicts(path_to_generated_output)
    # expected = read_csv_as_dicts(path_to_template_output)

    # assert actual == expected

def test_run_predict_structural_variants_unmasked():

    command = ["seabreeze", "run", "predict_structural_variants", "--dir","test/seabreeze_run_test","--assembly", "test/seabreeze_run_test/test_genomes/Ara+1_10000gen_4530A_truncated.fasta", "--ancestor", "test/seabreeze_run_test/test_genomes/Anc-_0gen_REL606_truncated.fasta", "--ori","AGGATGCTTTACCCAATATCAGCGAT","--dif","CTGGTCGGCAGAATGAGCAATCGCCA", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    # path_to_template_output="test/seabreeze_run_test/expected/Ara+1_10000gen_4530A_truncated_clean_unmasked.syri.out"
    # path_to_generated_output="test/seabreeze_run_test/07_syri_output/unmasked/Ara+1_10000gen_4530A_truncated/Ara+1_10000gen_4530A_truncated_clean.syri.out"

    # actual = read_csv_as_dicts(path_to_generated_output)
    # expected = read_csv_as_dicts(path_to_template_output)

    # assert actual == expected

def test_run_predict_structural_variants_masked():

    command = ["seabreeze", "run", "predict_structural_variants", "--dir","test/seabreeze_run_test","--masked","--assembly", "test/seabreeze_run_test/test_genomes/Ara+1_10000gen_4530A_truncated.fasta", "--ancestor", "test/seabreeze_run_test/test_genomes/Anc-_0gen_REL606_truncated.fasta", "--ori","AGGATGCTTTACCCAATATCAGCGAT","--dif","CTGGTCGGCAGAATGAGCAATCGCCA", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    # path_to_template_output="test/seabreeze_run_test/expected/Ara+1_10000gen_4530A_truncated_clean_masked.syri.out"
    # path_to_generated_output="test/seabreeze_run_test/07_syri_output/unmasked/Ara+1_10000gen_4530A_truncated/Ara+1_10000gen_4530A_truncated_clean.syri.out"

    # actual = read_csv_as_dicts(path_to_generated_output)
    # expected = read_csv_as_dicts(path_to_template_output)

    # assert actual == expected

def test_run_predict_replichore_balance_unmasked():

    command = ["seabreeze", "run", "predict_replichore_balance", "--dir","test/seabreeze_run_test","--assembly", "test/seabreeze_run_test/test_genomes/Ara+1_10000gen_4530A_truncated.fasta", "--ancestor", "test/seabreeze_run_test/test_genomes/Anc-_0gen_REL606_truncated.fasta", "--ori","AGGATGCTTTACCCAATATCAGCGAT","--dif","CTGGTCGGCAGAATGAGCAATCGCCA", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    # path_to_template_output="test/seabreeze_run_test/expected/inversion_replichores_unmasked.csv"
    # path_to_generated_output="test/seabreeze_run_test/11_annotated_boundaries/unmasked/inversion_replichores.csv"

    # actual = read_csv_as_dicts(path_to_generated_output)
    # expected = read_csv_as_dicts(path_to_template_output)

    # assert actual == expected

def test_run_predict_replichore_balance_masked():

    command = ["seabreeze", "run", "predict_replichore_balance", "--dir","test/seabreeze_run_test","--masked","--assembly", "test/seabreeze_run_test/test_genomes/Ara+1_10000gen_4530A_truncated.fasta", "--ancestor", "test/seabreeze_run_test/test_genomes/Anc-_0gen_REL606_truncated.fasta", "--ori","AGGATGCTTTACCCAATATCAGCGAT","--dif","CTGGTCGGCAGAATGAGCAATCGCCA", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    # path_to_template_output="test/seabreeze_run_test/expected/inversion_replichores_masked.csv"
    # path_to_generated_output="test/seabreeze_run_test/11_annotated_boundaries/masked/inversion_replichores.csv"


    # actual = read_csv_as_dicts(path_to_generated_output)
    # expected = read_csv_as_dicts(path_to_template_output)

    # assert actual == expected

def test_run_predict_SV_mechanism_unmasked():

    command = ["seabreeze", "run", "predict_SV_mechanism","--dir","test/seabreeze_run_test","--assembly", "test/seabreeze_run_test/test_genomes/Ara+1_10000gen_4530A_truncated.fasta", "--ancestor", "test/seabreeze_run_test/test_genomes/Anc-_0gen_REL606_truncated.fasta", "--ori","AGGATGCTTTACCCAATATCAGCGAT","--dif","CTGGTCGGCAGAATGAGCAATCGCCA", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"


def test_run_predict_SV_mechanism_masked():

    command = ["seabreeze", "run", "predict_SV_mechanism","--dir","test/seabreeze_run_test","--masked","--assembly", "test/seabreeze_run_test/test_genomes/Ara+1_10000gen_4530A_truncated.fasta", "--ancestor", "test/seabreeze_run_test/test_genomes/Anc-_0gen_REL606_truncated.fasta", "--ori","AGGATGCTTTACCCAATATCAGCGAT","--dif","CTGGTCGGCAGAATGAGCAATCGCCA", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

def test_run_annotate_SV_regions_unmasked():

    command = ["seabreeze", "run", "annotate_SV_regions","--dir","test/seabreeze_run_test","--assembly", "test/seabreeze_run_test/test_genomes/Ara+1_10000gen_4530A_truncated.fasta", "--ancestor", "test/seabreeze_run_test/test_genomes/Anc-_0gen_REL606_truncated.fasta", "--ori","AGGATGCTTTACCCAATATCAGCGAT","--dif","CTGGTCGGCAGAATGAGCAATCGCCA", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

def test_run_annotate_SV_regions_unmasked():

    command = ["seabreeze", "run", "annotate_SV_regions","--dir","test/seabreeze_run_test","--masked","--assembly", "test/seabreeze_run_test/test_genomes/Ara+1_10000gen_4530A_truncated.fasta", "--ancestor", "test/seabreeze_run_test/test_genomes/Anc-_0gen_REL606_truncated.fasta", "--ori","AGGATGCTTTACCCAATATCAGCGAT","--dif","CTGGTCGGCAGAATGAGCAATCGCCA", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

def test_remove_test_data():

    """ delete all test data.. """
    files_to_remove = ["test/seabreeze_run_test/03_reindex_genomes","test/seabreeze_run_test/04_rename_genome","test/seabreeze_run_test/05_isescan_tables","test/seabreeze_run_test/06_nucmer_alignment","test/seabreeze_run_test/07_syri_output","test/seabreeze_run_test/08_reindex_genome_oric","test/seabreeze_run_test/09_annotated_genomes","test/seabreeze_run_test/11_annotated_boundaries","test/seabreeze_run_test/12_genome_diff_tables"]

    for file in files_to_remove:
        shutil.rmtree(file)
