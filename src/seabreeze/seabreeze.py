"""
seabreeze.py, a wrapper script around a snakemake pipeline to analyse structural variants in bacterial genomes
Author: Ira Zibbu
Version: 0.0.3
"""

""" imports """
import click
import seabreeze
import importlib.resources
import subprocess
import os
import logging
import coloredlogs
from seabreeze import __version__
import shutil
import csv
import pandas as pd

# set up logger
coloredlogs.install(level='DEBUG')
logger = logging.getLogger('seabreeze')


def get_snakefile():

    """
    Fetch path to the seabreeze directory where the snakefile
    args: None
    returns: Path to the snakefile
    """

    seabreeze_package_path = importlib.resources.files(seabreeze)
    print(f"seabreeze package path is {seabreeze_package_path}")
    print(os.listdir(seabreeze_package_path))
    snakefile_path = os.path.join(seabreeze_package_path, "Snakefile")
    return snakefile_path

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)

def cli():

    """
    seabreeze is a tool for comprehensively analyzing genetic variation among bacterial genomes caused by structural mutations.
    To get help, run `seabreeze --help`

    Please see documentation at https://barrick.github.io/seabreeze
    Report bugs, errors and request features at https://github.com/barricklab/seabreeze

    seabreeze was developed at the Barrick Lab at UT Austin
    """

@cli.command(
    "batch",
    context_settings=dict(ignore_unknown_options=True),
    short_help="Run seabreeze in batch mode to process multiple samples at once",
)
@click.argument(
    "workflow",
    type=click.Choice(
        [
            "analyse_genome_sizes",
            "predict_IS_elements",
            "predict_structural_variants",
            "predict_replichore_balance",
            "predict_SV_mechanism",
            "annotate_SV_regions",
            "run_all"
        ]
    ),
)
@click.option(
    "-d",
    "--dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Set working directory. Default is current working directory",
    default=".",
)

@click.option(
    "-b",
    "--data",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Path to data.csv file. Default is data.csv in current working dir",
    default="./data.csv",
)

@click.option(
    "--oridif",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Path to csv file with sequences for the origin and terminus. Default = ori_dif_coords.csv",
    default="./ori_dif_coords.csv",
)

@click.option(
    "--threads",
    type=int,
    default=4,
    show_default=True,
    help="Resources to be used.",
)

@click.option("--masked", is_flag=True, help="Mask insertion sequences")

@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)

def run_smk_batch_workflow(workflow,dir,data,oridif,masked,threads,snakemake_args):

    """ Determine snakemake command to run based on supplied args """

    logger.debug(f"==================================")
    logger.debug(f"seabrezee version: {__version__}")
    logger.debug(f"==================================")

    # make the working directory if it does not already exist
    if not os.path.exists(dir):
        logger.debug(f"{dir} does not exist. Creating the directory..")
        try:
            os.makedirs(dir)
        except OSError as e:
            logger.critical(f"Unable to create {dir}")
            logger.critical(e)

    # verify that data.csv exists
    if not os.path.exists(data):
        logger.critical(f"{data} does not exist. Aborting")
        exit(1)


    df = pd.read_csv(data)

    # Check if the csv has exactly two columns
    if df.shape[1] != 2:
        logger.critical("data.csv must have exactly two columns.")

    # Check if the columns are named correctly
    expected_columns = ['assembly', 'ancestor']
    if list(df.columns) != expected_columns:
        logger.critical(f"The data.csv columns must be named {expected_columns}. Found: {list(df.columns)}")

    # Check if there are any duplicate entries for assemblies
    if df['assembly'].duplicated().any():
        logger.critical("There are duplicate entries in the assembly column.")

    list_of_ancestors = df['ancestor'].tolist()
    list_of_assemblies = df['assembly'].tolist()

    for ancestor in list_of_ancestors:
        if not(ancestor in list_of_assemblies):
            logger.critical(f"Ancestor genome {ancestor} needs a self-to-self comparison. Please append this line to data.csv: {ancestor},{ancestor}")


    # ori and dif sequences are required to predict replichore balance:
    if workflow == "predict_replichore_balance":
        if not os.path.exists(oridif):
            logger.critical(f"{oridif} does not exist. Aborting")
            exit(1)

    os.chdir(dir)
    logger.debug(f"Switch working directory to {os.getcwd()}")
    logger.debug(f"seabrezee version: {__version__}")

    masked_workflows=["predict_structural_variants","predict_replichore_balance","predict_SV_mechanism","annotate_SV_regions","run_all"]
    if workflow in masked_workflows and masked:
        workflow=f"{workflow}_masked"

    elif workflow in masked_workflows and not(masked):
        workflow=f"{workflow}_unmasked"

    seabreeze_package_path = importlib.resources.files(seabreeze)

    config_options=f"data={data} oridif={oridif}"
    cmd=()

    cmd = (
        "snakemake --snakefile {snakefile} --use-conda --cores {cores} {workflow} --config {config}"
        " {args} "
    ).format(
        snakefile=get_snakefile(),
        cores=threads,
        workflow=workflow,
        args=" ".join(snakemake_args),
        config=config_options
    )
    logger.debug("Executing: %s" % cmd)

    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e)
        exit(1)

# ---- SINGLE RUN MODE ----

@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True),
    short_help="Run seabreeze for a single sample pair",
)

@click.argument(
    "workflow",
    type=click.Choice(
        [
            "analyse_genome_sizes",
            "predict_IS_elements",
            "predict_structural_variants",
            "predict_replichore_balance",
            "predict_SV_mechanism",
            "annotate_SV_regions",
            "run_all"
        ]
    ),
)
@click.option(
    "-d",
    "--dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Set working directory. Default is current working directory",
    default=".",
)

@click.option(
    "--ancestor",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Path to FASTA file of the ancestor or reference sequence",
)

@click.option(
    "--assembly",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Path to FASTA file of the assembly or query sequence",
)

@click.option(
    "--threads",
    type=int,
    default=4,
    show_default=True,
    help="Resources to be used.",
)

@click.option(
    "--ori",
    required=False,
    help="Sequence of the origin. Must use IUPAC bases, and be unique in both the ancestor and assembly",
)

@click.option(
    "--dif",
    required=False,
    help="Sequence of the origin. Must use standard IUPAC bases, and be unique in both the ancestor and assembly",
)

@click.option("--masked", is_flag=True, help="Mask insertion sequences")

@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)

def run_smk_single_workflow(workflow,dir,assembly,ancestor,ori,dif,masked,threads,snakemake_args):

    """ Determine snakemake command to run based on supplied args """

    logger.debug(f"==================================")
    logger.debug(f"seabrezee version: {__version__}")
    logger.debug(f"==================================")


    # check if working dir exists. If not, create.

    if not os.path.exists(dir):
        logger.debug(f"{dir} does not exist. Creating the directory..")
        try:
            os.makedirs(dir)
        except OSError as e:
            logger.critical(f"Unable to create {dir}")
            logger.critical(e)
            exit(1)

    # check that supplied fasta files exist
    for file in [ancestor,assembly]:
        if not os.path.exists(file):
            logger.critical(f"{file} does not exist. Aborting")
            exit(1)

    dir_for_fasta_files = os.path.join(dir, "02_genomes")

    # make the dir for the fasta files and copy them over

    if os.path.exists(dir_for_fasta_files):
        logger.debug(f"{dir_for_fasta_files} exists. No requirement to make file")

    if not os.path.exists(dir_for_fasta_files):
        try:
            os.makedirs(dir_for_fasta_files)
        except OSError as e:
            logger.critical(f"Unable to create {dir_for_fasta_files}")
            logger.critical(e)
            exit(1)

    try:
        for file_path in [ancestor, assembly]:
            shutil.copy(file_path, dir_for_fasta_files)
            logger.debug(f"Copied {file_path} to {dir_for_fasta_files}")
    except Exception as e:
        logger.critical(f"Error copying files: {e}")
        exit(1)

    # define the data.csv from scratch

    os.chdir(dir)
    logger.debug(f"Switch working directory to {os.getcwd()}")

    ancestor_name = os.path.splitext(os.path.basename(ancestor))[0]
    assembly_name = os.path.splitext(os.path.basename(assembly))[0]

    data_csv_table = [
        {'assembly': assembly_name, 'ancestor': ancestor_name},
        {'assembly': ancestor_name, 'ancestor': ancestor_name}
    ]

    csv_file = "data.csv"
    try:
        with open(csv_file, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=['ancestor', 'assembly'])
            writer.writeheader()
            writer.writerows(data_csv_table)
        logger.debug(f"data.csv written successfully to {csv_file}")
    except Exception as e:
        logger.critical(f"Error writing CSV: {e}")

    # ori and dif sequences are required to predict replichore balance:
    if workflow == "predict_replichore_balance" or workflow == "run_all":
        if not ori or not dif:
            logger.critical("Both --ori and --dif options must be supplied for the predict_replichore_balance or run_all workflow.")
            exit(1)

        # write out an ori_dif_sequences.csv file

    ori_dif_table = [
        {'ancestor': ancestor_name, 'ori': ori, 'dif': dif}]

    csv_file = "ori_dif_sequences.csv"
    try:
        with open(csv_file, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=['ancestor', 'ori', 'dif'])
            writer.writeheader()
            writer.writerows(ori_dif_table)
        logger.debug(f"ori_dif_sequences.csv written successfully to {csv_file}")
    except Exception as e:
        logger.critical(f"Error writing CSV: {e}")
        exit(1)

    masked_workflows=["predict_structural_variants","predict_replichore_balance","predict_SV_mechanism","annotate_SV_regions","run_all"]
    if workflow in masked_workflows and masked:
        workflow=f"{workflow}_masked"

    elif workflow in masked_workflows and not(masked):
        workflow=f"{workflow}_unmasked"

    seabreeze_package_path = importlib.resources.files(seabreeze)

    config_options=f"data=data.csv oridif=ori_dif_sequences.csv"
    cmd=()

    cmd = (
        "snakemake --snakefile {snakefile} --use-conda --cores {cores} {workflow} --config {config}"
        " {args} "
    ).format(
        snakefile=get_snakefile(),
        cores=threads,
        workflow=workflow,
        args=" ".join(snakemake_args),
        config=config_options
    )
    logger.debug("Executing: %s" % cmd)

    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logger.critical(e)
        exit(1)

def main():
    cli()

if __name__ == "__main__":
    main()
