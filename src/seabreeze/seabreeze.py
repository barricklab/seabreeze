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
import logging
from seabreeze import __version__

# set up logger
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('seabreeze')


def get_snakefile():

    """
    Fetch path to the seabreeze directory where the snakefile
    args: None
    returns: Path to the snakefile
    """

    seabreeze_package_path = importlib.resources.files(seabreeze)
    snakefile_path = os.path.join(seabreeze_package_path, "workflow/Snakefile")
    return snakefile_path



@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)

def cli():
    """
    This is information about seabreeze
    """

@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True),
    short_help="Run seabreeze",
)
@click.argument(
    "workflow",
    type=click.Choice(
        [
            "analyse_genome_sizes",
            "predict_IS_elements",
            "predict_structural_variants",
            "copy_fasta_files",
            "all"
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
    help="Resources to be used. Default = 4. `all` is a valid option",
)

@click.option("--masked", is_flag=True, help="Mask insertion sequences")
@click.option("--batch", is_flag=True, help="Enable batch mode. Requires --data argument to be supplied")

@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)

def run_smk_workflow(workflow,dir,data,oridif, ancestor,assembly,masked,batch,threads,snakemake_args):

    os.chdir(dir)
    logger.info(f"seabrezee version: {__version__}")

    masked_workflows=["predict_structural_variants","predict_replichore_balance","predict_SV_mechanism","annotate_SV_regions"]
    if workflow in masked_workflows and masked:
        workflow=f"{workflow}_masked"

    seabreeze_package_path = importlib.resources.files(seabreeze)

    config_options=f"data={data}"
    cmd=()

    if batch:
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
        # removes the traceback
        logger.critical(e)
        exit(1)

def main():
    cli()

if __name__ == "__main__":
    main()
