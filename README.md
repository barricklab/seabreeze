# Seabreeze
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

3. 
### Misc information
Where does seabreeze get its name from? **S**tructural **v**ariant for **b**acterial **rese**quencing abbreviates as SVbrese, which sounded like seabreeze. 