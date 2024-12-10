# Preparing input

Assuming you have [completed installation](installation.md), you can now prepare input files for analysis. At its core _seabreeze_ is designed for bacterial resequencing experiments, and therefore predicts mutations between reference-query pairs, where the reference sequence is referred to as an _ancestor_ and the query sequence being compared to is called the _assembly_. This nomenclature is use throughout the documentation and the software.
_seabreeze_ requires FASTA files for all of the bacterial genomes being analysed. These FASTA files should only contain a single contig for the full-length genome. It should not contain additional entries for plasmids. _seabreeze_ works best with genome assemblies generated from long-read sequencing as these can capture structural variants that occur between repeat sequences. It also does not work well for comparing genomes with high sequence/ phylogenetic divergence (for example between different bacterial species).

All of the input/output of the analysis is stored in `data/`. The pairwise comparisons to be performed are specified in a csv file called `data.csv` in the `data/` directory. Populate the csv file with the ancestor-assembly pairs you wish to analyse. For example:

| assembly | ancestor |
| -------|---------|
| genome1 | genome2 |
| genome3 | genome4 |

Place the corresponding FASTA files for these genomes in `02_genomes`. _seabreeze_ expects that these file names to match, and contain the extension `.fasta`. Therefore, for this example, the `data/` directory should look like this:
```
seabreeze
|
|---data/
|   |
|   |---data.csv
|   |---02_genomes/
|   |   |
|   |   |---genome1.fasta
|   |   |---genome2.fasta
|   |   |---genome3.fasta
|   |   |---genome4.fasta
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

Some commands ahead may require additional input, and will be specified in respective sections.

# Commands

All of the _seabreeze_ commands take the following format and should be executed from the _seabreeze_ root directory:
```
snakemake --use-conda --cores n <command>
```
Where n is the number of cores to be used to execute the command. Following sections describe commands that can be used. Additional flags can be added to this command and are detailed in 'Additional Options`.

## Analyse genome sizes

This command compares the sizes of the assemblies to their specified ancestors in a pairwise comparisons. 

```
snakemake --use-conda --cores n analyse_genome_sizes
```
The output file generated is `data/04_rename_genome/genome_size_stats.csv`. Please see [Output](output.md) for more information about the fields in this table.


