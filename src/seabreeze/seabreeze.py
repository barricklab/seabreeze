"""
seabreeze.py, a wrapper script around a snakemake pipeline to analyse structural variants in bacterial genomes
Author: Ira Zibbu
Version: 0.0.1
"""

""" imports """
import click
import seabreeze
import importlib.resources
import subprocess
import os
from seabreeze import __version__

def get_package_dir():

    """
    Fetch path to the seabreeze directory where the snakemake rules live
    args: None
    returns: Path to the directory with the snakemake files
    """

    seabreeze_package_path = importlib.resources.files(seabreeze)
    rules_path = os.path.join(seabreeze_package_path, "workflow", "rules")

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)

def cli():
    """
    This is information about seabreeze
    """

@cli.command()
def greet():
    """Say hello"""
    click.echo("Hello!")

def main():
    cli()

if __name__ == "__main__":
    main()
