# Installation

## Requirements
- A UNIX-based operating system (like Linux or MacOS) with bash is recommended. This software has not been tested with the Windows shell. 
- _seabreeze_ requires conda >= 23.3.1 to install python packages and manage environments. We recommend [Miniforge](https://github.com/conda-forge/miniforge), a minimal conda installer that is similar to [Miniconda](https://docs.anaconda.com/miniconda/). 

## Recommended: Install via bioconda

We recommend creating a dedicated conda environment for seabreeze, to avoid package conflicts. 

```
conda env create --name seabreeze
```

Once built, activate the environment:

```
conda activate seabreeze
```

Finally, install the package from bioconda.

```
conda install seabreeze-genomics -c bioconda
```

In case this does not work, try specifying:

```
conda install seabreeze-genomics -c conda-forge -c bioconda
```


## Install from GitHub


Download the latest release of _seabreeze_ from the [Github Repository](https://github.com/barr>

_seabreeze_ only requires an intial conda environment to be set up. Navigate to the _seabreeze_ diectory and run the following command:

```
 conda env create --name seabreeze --file environment.yml
```

Activate the conda environment before running the commands.

```
conda activate seabreeze
```

Finally,

```
pip install .
```

You're all set! See [Usage](usage.md) to get started.

## Optional
By default, _seabreeze_ creates conda environments and installs packages in the _seabreeze_  directory. However, if you plan to use _seabreeze_ multiple times, or with different directories for different samples, then it can be helpful to have a common location on your system to store environments. This can save disk space and time by avoiding repeated installation of environments. 

Create a folder for Snakemake to store conda environments. We recommend a directory in $HOME, but it could be anywhere.
```
mkdir $HOME/snakemake_conda_envs
```

Append this line to your `~/.bashrc`, `~/.zshrc` or as applicable to the shell you are using, using the path to the directory created in the previous step.
```
export SNAKEMAKE_CONDA_PREFIX=$HOME/snakemake_conda_envs
```
   


