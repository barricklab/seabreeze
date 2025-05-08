"""
Master test to see if the whole seabreeze pipeline runs
"""

import subprocess
import pytest
import csv

def test_analyse_genome_sizes():

    command = ["seabreeze", "run", "--batch", "analyse_genome_sizes", "--data", "test/data/seabreeze_test/data.csv", "--dir", "test/data/seabreeze_test/"]
    result = subprocess.run(command)
    assert result.returncode == 0, f"Command failed with error: {result.stderr}"

    path_to_template_output="test/data/seabreeze_test/expected/genome_size_stats.csv"
    path_to_generated_output="test/data/seabreeze_test/04_rename_genome/genome_size_stats.csv"

    actual = read_csv_as_list_of_dicts(path_to_generated_output)
    expected = read_csv_as_list_of_dicts(path_to_template_output)

    assert actual == expected
