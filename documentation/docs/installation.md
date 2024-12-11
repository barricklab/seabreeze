# Installation

## Requirements
- A UNIX-based operating system (like Linux or MacOS) with bash is recommended. This software has not been tested with the Windows shell. 
- _seabreeze_ requires conda >= 24.3.0 to install python packages and manage environments. We recommend [Miniforge](https://github.com/conda-forge/miniforge), a minimal conda installer that is similar to [Miniconda](https://docs.anaconda.com/miniconda/). 

## Installation
_seabreeze_ is a Snakemake workflow with associated Python scripts. Download the latest release of _seabreeze_ from the [Github Repository](https://github.com/barricklab/seabreeze/releases).

## Set up
_seabreeze_ only requires an intial conda environment to be set up. Navigate to the _seabreeze_ diectory and run the following command:
```
 conda env create --name seabreeze --file bin/workflow/envs/seabreeze.yml
```

Activate the conda environment before running the commands

```
conda activate seabreeze
```

You're all set! See [Usage](usage.md) to get started.

## Optional
By default, _seabreeze_ creates conda environments and installs packages in the _seabreeze_ working directory. However, if you plan to use _seabreeze_ multiple times, or with different directories for different samples, then it can be helpful to have a common location on your system to store environments. This can save disk space and time. 

Create a folder for Snakemake to store conda environments. We recommend a directory in $HOME, but it could be anywhere.
```
mkdir $HOME/snakemake_conda_envs
```

Append this line to your `~/.bashrc`, `~/.zshrc` or as applicable to the shell you are using:
```
export SNAKEMAKE_CONDA_PREFIX=$HOME/snakemake_conda_envs
```
   


