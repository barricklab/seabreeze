"""
Master test to see if the whole seabreeze pipeline runs
"""

import subprocess
import pytest
import csv

def read_csv_as_dicts(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))

def test_analyse_genome_sizes():

    command = ["seabreeze", "batch", "analyse_genome_sizes", "--data", "test/seabreeze_test/data.csv", "--dir", "test/seabreeze_test/"]
    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    path_to_template_output="test/seabreeze_test/expected/genome_size_stats.csv"
    path_to_generated_output="test/seabreeze_test/04_rename_genome/genome_size_stats.csv"

    actual = read_csv_as_dicts(path_to_generated_output)
    expected = read_csv_as_dicts(path_to_template_output)

    assert actual == expected

def test_predict_IS_elements():

    command = ["seabreeze", "batch", "--dir","test/seabreeze_test/","--data", "test/seabreeze_test/data.csv", "predict_IS_elements"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    path_to_template_output="test/seabreeze_test/expected/Anc-_0gen_REL606_truncated.csv"
    path_to_generated_output="test/seabreeze_test/05_isescan_tables/Anc-_0gen_REL606_truncated.csv"

    actual = read_csv_as_dicts(path_to_generated_output)
    expected = read_csv_as_dicts(path_to_template_output)

    assert actual == expected

def test_predict_structural_variants_unmasked():

    command = ["seabreeze", "batch", "--dir","test/seabreeze_test/","--data", "test/seabreeze_test/data.csv", "predict_structural_variants"]

    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    path_to_template_output="test/seabreeze_test/expected/Ara+1_10000gen_4530A_truncated_clean.syri.out"
    path_to_generated_output="test/seabreeze_test/07_syri_output/unmasked/Ara+1_10000gen_4530A_truncated/Ara+1_10000gen_4530A_truncated_clean.syri.out"

    actual = read_csv_as_dicts(path_to_generated_output)
    expected = read_csv_as_dicts(path_to_template_output)

    assert actual == expected
