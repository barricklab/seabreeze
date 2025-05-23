# Preparing input

Assuming you have [completed installation](installation.md), you can now prepare input files for analysis. At its core, _seabreeze_ is designed for bacterial resequencing experiments, and therefore predicts mutations between reference-query pairs, where the reference sequence is referred to as the _ancestor_ and the query sequence being compared to is called the _assembly_. This nomenclature is use throughout the documentation and the software.
_seabreeze_ requires FASTA files for all of the bacterial genomes being analysed. These FASTA files should only contain a single contig for the full-length genome. It should not contain additional entries for plasmids. _seabreeze_ works best with genome assemblies generated from long-read sequencing as these can capture structural variants that occur between repeat sequences. It also does not work well for comparing genomes with high sequence/ phylogenetic divergence (for example between different bacterial species).

_seabreeze_ has two modes: 
- `run`:  Analyzing one ancestor-assembly pair
- `batch`: Analyzing multiple samples at once

Both modes run the same workflows, but require input to be supplied differently. 
## `run`: Analyzing one ancestor-assembly pair

```
Usage: seabreeze run [OPTIONS] [WORKFLOW] [SNAKEMAKE_ARGS]...

Options:
  -d, --dir PATH     Set working directory. Default is current working
                     directory
  --ancestor PATH    Path to FASTA file of the ancestor or reference sequence
  --assembly PATH    Path to FASTA file of the assembly or query sequence
  --threads INTEGER  Resources to be used.  [default: 4]
  --ori TEXT         Sequence of the origin. Must use IUPAC bases, and be
		             unique in both the ancestor and assembly.  Required for                          predict_replichore_balance or run_all workflows
  --dif TEXT         Sequence of the terminus. Must use standard IUPAC bases,
                     and be unique in both the ancestor and assembly. Required                        for predict_replichore_balance or run_all workflows
  --masked           Mask insertion sequences
  -h, --help         Show this message and exit.

```

- All input is directly supplied via the command line. No set up required.
## `batch`: Analyzing multiple samples at once

```
Usage: seabreeze batch [OPTIONS] [WORKFLOW] [SNAKEMAKE_ARGS]...

  Determine snakemake command to run based on supplied args

Options:
  -d, --dir PATH     Set working directory. Default is current working
                     directory
  -b, --data PATH    Path to data.csv file. Default is data.csv in current
                     working dir
  --oridif PATH      Path to csv file with sequences for the origin and
                     terminus. Default is ori_dif_coords.csv in the current dir.
                     Required for predict_replichore_balance or run_all workflows
  --threads INTEGER  Resources to be used.  [default: 4]
  --masked           Mask insertion sequences
  -h, --help         Show this message and exit.

```

`batch` mode requires the following set up:

- `data.csv`: The pairwise comparisons to be performed are specified in a csv file called `data.csv` , whose path is supplied via the `--data` option. This file has only two columns `ancestor` and `assembly`, in a standard csv format. Populate the csv file with the ancestor-assembly pairs you wish to analyse. For example:

| assembly | ancestor |
| -------- | -------- |
| genome1  | genome2  |
| genome3  | genome4  |

- `02_genomes`: A subdirectory of the specified working directory in the `--dir` option. Place the corresponding FASTA files for these genomes in `02_genomes`. _seabreeze_ expects that these file names to match, and contain the extension `.fasta`. For example the corresponding working directory should look like this

```
<my_working_dir>
|
|---02_genomes/
|   |
|   |---genome1.fasta
|   |---genome2.fasta
|   |---genome3.fasta
|   |---genome4.fasta
```

Note that many assemblies can use the same ancestor, but the same assembly cannot use many ancestors:
For example, this is fine:
 
| assembly | ancestor |
|-------|-------|
| genome1 | genome2 |
| genome3 | genome2 |

This will cause a fatal error:

| assembly | ancestor |
|--------|-------|
| genome1 | genome2 |
| genome1 | genome3 |

It is also a good idea to set up a comparison of the ancestor genomes to themselves. This can be a quality control step, as we do not expect the pipeline to predict any mutations in a genome relative to itself. The following table demonstrates an implementation of this:

| assembly | ancestor |
| -------- | -------- |
| genome1  | genome2  |
| genome3  | genome4  |
| genome2  | genome2  |
| genome4  | genome4  |

---

## Valid workflows

- `analyse_genome_sizes`
- `predict_IS_elements`
- `predict_structural_variants`
- `predict_replichore_balance`
- `predict_SV_mechanism`
- `annotate_SV_regions`
- `run_all`

For detailed information on each workflow, please see [the workflows page](workflows.md)

---
## Options

These options can be 

- `--masked`: Mask insertion sequences with NNNs. Masking repeats can help detect structural mutations.

---
## Additional Snakemake argument

Any valid Snakemake argument can be passed at the end of the command, but the following are the most commonly used options

- `-np`: This flag performs a 'dry-run' of the command i.e. it checks to see that all of the input files necessary for that command are present, or can be generated by _seabreeze_. It can helpful to run this before a computationally expensive step such as predicting IS elements, or predicting SVs.

If a command was terminated before it was completed, Snakemake may prompt you to use these flags when you run the same command again:

- `--unlock`:  Snakemake will "lock" a directory if all of the output files were not completed before the processes were terminated. Use this flag to unlock the directory.
- `--rerun-incomplete`:  To allow Snakemake to remake any output files that were not completed in the prior run.

- `--conda-frontend conda`: Some users have reported a `CreatedCondaEnvironmentException` that arises when trying to run any of the Snakemake workflows (see below). Appending this flag should address this issue. 

```
CreateCondaEnvironmentException:
The 'mamba' command is not available in the shell /usr/bin/bash that will be used by Snakemake. You have to ensure that it is in your PATH, e.g., first activating the conda base environment with `conda activate base`.The mamba package manager (https://github.com/mamba-org/mamba) is a fast and robust conda replacement. It is the recommended way of using Snakemake's conda integration. It can be installed with `conda install -n base -c conda-forge mamba`. If you still prefer to use conda, you can enforce that by setting `--conda-frontend conda`.
```

