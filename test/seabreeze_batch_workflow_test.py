"""
Master test to see if the seabreeze workflow pipelines run
"""

import subprocess
import pytest
import csv
import os
import shutil

def read_csv_as_dicts(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))

def test_batch_analyse_genome_sizes():

    command = ["seabreeze", "batch", "analyse_genome_sizes", "--data", "test/seabreeze_batch_test/data.csv", "--dir", "test/seabreeze_batch_test/", "--conda-frontend", "conda"]
    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    path_to_template_output="test/seabreeze_batch_test/expected/genome_size_stats.csv"
    path_to_generated_output="test/seabreeze_batch_test/04_rename_genome/genome_size_stats.csv"

    # actual = read_csv_as_dicts(path_to_generated_output)
    # expected = read_csv_as_dicts(path_to_template_output)

    # assert actual == expected

def test_batch_predict_IS_elements():

    command = ["seabreeze", "batch", "--dir","test/seabreeze_batch_test/","--data", "test/seabreeze_batch_test/data.csv", "predict_IS_elements", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    path_to_template_output="test/seabreeze_batch_test/expected/Anc-_0gen_REL606_truncated.csv"
    path_to_generated_output="test/seabreeze_batch_test/05_isescan_tables/Anc-_0gen_REL606_truncated.csv"

    # actual = read_csv_as_dicts(path_to_generated_output)
    # expected = read_csv_as_dicts(path_to_template_output)

    # assert actual == expected

def test_batch_predict_structural_variants_unmasked():

    command = ["seabreeze", "batch", "--dir","test/seabreeze_batch_test/","--data", "test/seabreeze_batch_test/data.csv", "predict_structural_variants", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    # path_to_template_output="test/seabreeze_batch_test/expected/Ara+1_10000gen_4530A_truncated_clean_unmasked.syri.out"
    # path_to_generated_output="test/seabreeze_batch_test/07_syri_output/unmasked/Ara+1_10000gen_4530A_truncated/Ara+1_10000gen_4530A_truncated_clean.syri.out"

    # actual = read_csv_as_dicts(path_to_generated_output)
    # expected = read_csv_as_dicts(path_to_template_output)

    # assert actual == expected

def test_batch_predict_structural_variants_masked():

    command = ["seabreeze", "batch", "--dir","test/seabreeze_batch_test/","--data", "test/seabreeze_batch_test/data.csv", "--masked", "predict_structural_variants", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    # path_to_template_output="test/seabreeze_batch_test/expected/Ara+1_10000gen_4530A_truncated_clean_masked.syri.out"
    # path_to_generated_output="test/seabreeze_batch_test/07_syri_output/unmasked/Ara+1_10000gen_4530A_truncated/Ara+1_10000gen_4530A_truncated_clean.syri.out"

    # actual = read_csv_as_dicts(path_to_generated_output)
    # expected = read_csv_as_dicts(path_to_template_output)

    # assert actual == expected

def test_batch_predict_replichore_balance_unmasked():

    command = ["seabreeze", "batch", "--dir","test/seabreeze_batch_test/","--data", "test/seabreeze_batch_test/data.csv", "--oridif", "test/seabreeze_batch_test/ori_dif_sequences.csv", "predict_replichore_balance", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    # path_to_template_output="test/seabreeze_batch_test/expected/inversion_replichores_unmasked.csv"
    # path_to_generated_output="test/seabreeze_batch_test/11_annotated_boundaries/unmasked/inversion_replichores.csv"

    # actual = read_csv_as_dicts(path_to_generated_output)
    # expected = read_csv_as_dicts(path_to_template_output)

    # assert actual == expected

def test_batch_predict_replichore_balance_unmasked():

    command = ["seabreeze", "batch", "--dir","test/seabreeze_batch_test/","--data", "test/seabreeze_batch_test/data.csv", "--masked", "--oridif", "test/seabreeze_batch_test/ori_dif_sequences.csv", "predict_replichore_balance", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    # path_to_template_output="test/seabreeze_batch_test/expected/inversion_replichores_masked.csv"
    # path_to_generated_output="test/seabreeze_batch_test/11_annotated_boundaries/masked/inversion_replichores.csv"


    # actual = read_csv_as_dicts(path_to_generated_output)
    # expected = read_csv_as_dicts(path_to_template_output)

    # assert actual == expected

def test_batch_predict_SV_mechanism_unmasked():

    command = ["seabreeze", "batch", "--dir","test/seabreeze_batch_test/","--data", "test/seabreeze_batch_test/data.csv", "--oridif", "test/seabreeze_batch_test/ori_dif_sequences.csv", "predict_SV_mechanism", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"


def test_batch_predict_SV_mechanism_masked():

    command = ["seabreeze", "batch", "--dir","test/seabreeze_batch_test/","--masked","--data", "test/seabreeze_batch_test/data.csv", "--oridif", "test/seabreeze_batch_test/ori_dif_sequences.csv", "predict_SV_mechanism", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

def test_batch_annotate_SV_regions_unmasked():

    command = ["seabreeze", "batch", "--dir","test/seabreeze_batch_test/","--data", "test/seabreeze_batch_test/data.csv", "--oridif", "test/seabreeze_batch_test/ori_dif_sequences.csv", "annotate_SV_regions", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

def test_batch_annotate_SV_regions_unmasked():

    command = ["seabreeze", "batch", "--dir","test/seabreeze_batch_test/","--data", "test/seabreeze_batch_test/data.csv", "--oridif", "test/seabreeze_batch_test/ori_dif_sequences.csv","--masked", "annotate_SV_regions", "--conda-frontend", "conda"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

def test_remove_test_data():

    """ delete all test data.. """
    files_to_remove = ["test/seabreeze_batch_test/03_reindex_genomes","test/seabreeze_batch_test/04_rename_genome","test/seabreeze_batch_test/05_isescan_tables","test/seabreeze_batch_test/06_nucmer_alignment","test/seabreeze_batch_test/07_syri_output","test/seabreeze_batch_test/08_reindex_genome_oric","test/seabreeze_batch_test/09_annotated_genomes","test/seabreeze_batch_test/11_annotated_boundaries","test/seabreeze_batch_test/12_genome_diff_tables"]

    for file in files_to_remove:
        shutil.rmtree(file)
