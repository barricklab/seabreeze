# Seabreeze
Author: Ira Zibbu
Last update to README: 2024-03-11

Seabreeze is an automated pipeline for structural variant analysis from bacterial resequencing experiments. It is currently under development as I migrate and clean the workflow first described in the [LTEE-SV](https://github.com/ira-zibbu/LTEE-SV) project repo. 

Seabreeze uses snakemake to neatly manage input/output file and conda environment dependencies.

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
mkdir $HOME/snakemake\_conda\_envs
```
Add this line to your .bashrc .zshrc or whatever shell you are currently using:
```
export SNAKEMAKE\_CONDA\_PREFIX=$HOME/snakemake\_conda\_envs
```

4. You're all set!
Navigate to the seabreeze directory and run snakemake with a few magic words:
```
conda activate seabreeze
snakemake --use-conda --cores 4 # feel free to change this but anything less than 2 cores will probably take an hour 
```

### To do:
Add testing data
 
### Misc information
Where does seabreeze get its name from? **S**tructural **v**ariant for **b**acterial **rese**quencing abbreviates as SVbrese, which sounded like seabreeze. 
