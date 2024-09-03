# Seabreeze
Author: Ira Zibbu

Last update to README: 2024-08-28

Seabreeze is an automated pipeline for structural variant analysis from bacterial resequencing experiments, and this repository holds the a version of seabreeze adapted to analyse structural variants

Seabreeze uses snakemake to neatly manage input/output file and conda environment dependencies.

### Updates
As of August 2024, this is the general purpose structural variant analysis tool [seabreeze](https://github.com/barricklab/seabreeze). The modified version of _seabreeze_ is now hosted separately on the my personal GitHub. Please use this repository's version for your own data, and not the modified version. 

# Quick start
To replicate analysis, follow these steps:

1. Clone this repository
```
git clone https://github.com/ira-zibbu/seabreeze.git
```


2. Create the seabreeze conda environment
```
 conda env create --name seabreeze --file bin/workflow/envs/seabreeze.yml
```


3. Create a folder for the conda environments

To allow Snakemake to keep reusing conda environments (i.e to prevent Snakemake from reinstalling environments again and again), do the following:
```
mkdir $HOME/snakemake_conda_envs
```
Add this line to your .bashrc .zshrc or whatever shell you are currently using:
```
export SNAKEMAKE_CONDA_PREFIX=$HOME/snakemake_conda_envs
```

4. Add test data

This is a small test data set to make sure everything runs okay
```
cp -r test data
cp data/data.csv ../data.csv
```


5. You're all set!
Navigate to the seabreeze directory and run snakemake with a few magic words:
```
conda activate seabreeze
snakemake --use-conda --cores 4 # feel free to change this but anything less than 2 cores will probably take an hour 
```

### Misc information
Where does seabreeze get its name from? **S**tructural **v**ariant for **b**acterial **rese**quencing abbreviates as SVbrese, which sounded like seabreeze. 
