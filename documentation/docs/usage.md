# Preparing input

Assuming you have [completed installation](installation.md), you can now prepare input files for analysis. At its core, _seabreeze_ is designed for bacterial resequencing experiments, and therefore predicts mutations between reference-query pairs, where the reference sequence is referred to as an _ancestor_ and the query sequence being compared to is called the _assembly_. This nomenclature is use throughout the documentation and the software.
_seabreeze_ requires FASTA files for all of the bacterial genomes being analysed. These FASTA files should only contain a single contig for the full-length genome. It should not contain additional entries for plasmids. _seabreeze_ works best with genome assemblies generated from long-read sequencing as these can capture structural variants that occur between repeat sequences. It also does not work well for comparing genomes with high sequence/ phylogenetic divergence (for example between different bacterial species).

All of the input/output of the analysis is stored in the `data/` directory in the _seabreeze_ folder. The pairwise comparisons to be performed are specified in a csv file called `data.csv` in the `data/` directory. This file has only two columns `ancestor` and `assembly`, in a standard csv format. Populate the csv file with the ancestor-assembly pairs you wish to analyse. For example:

| assembly | ancestor |
| -------- | -------- |
| genome1  | genome2  |
| genome3  | genome4  |

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

It is also a good idea to set up a comparison of the ancestor genomes to themselves. This can be a quality control step, as we do not expect the pipeline to predict any mutations in a genome relative to itself. The following table demonstrates an implementation of this:

| assembly | ancestor |
| -------- | -------- |
| genome1  | genome2  |
| genome3  | genome4  |
| genome2  | genome2  |
| genome4  | genome4  |

Some commands ahead may require additional input, and will be specified in respective sections.

---
# Commands

All of the _seabreeze_ commands take the following format and should be executed from the _seabreeze_ root directory:
```
snakemake --use-conda --cores n <command>
```
Where n is the number of cores to be used to execute the command. Following sections describe commands that can be used. Additional flags can be added to commands and are detailed in 'Additional Options'.

## Analyse genome sizes

This command compares the sizes of the assemblies to their specified ancestors in a pairwise manner. 

```
snakemake --use-conda --cores n analyse_genome_sizes
```
The output file generated is `data/04_rename_genome/genome_size_stats.csv`. Please see [output](output.md) for more information about the fields in this table, and the [tutorial](tutorial.md) for a sample.
## Find insertion sequences

This command uses [ISEscan](https://github.com/xiezhq/ISEScan/blob/master/README.md?plain=1) to annotate partial and full-length IS elements in the genomes. 
```
snakemake --use-conda --cores n predict_IS_elements
```

This command primarily generates a series of csv files`data/05_isescan_tables/<assembly>.csv`, with each csv file corresponding to an assembly in `data/data.csv`. ISEScan also generates additional output for each assembly, which is stored in their respective directories `data/05_iescan_tables/<assembly>/`.  For information about interpreting this csv file, we refer users to the official [ISEScan documentation](https://github.com/xiezhq/ISEScan/blob/master/README.md). 
## Predict SVs

This command aligns the assembly-ancestor pairs, predicts structural variants and generates synteny plots. 

```
snakemake --use-conda --cores n predict_structural_variants
```

- `data/07_syri_output/<assembly>/<assembly>_clean.syri.out` is a tsv file that describes the location and type of SVs in the assembly, relative to its ancestor. We refer the user to the official [SyRI documentation](https://schneebergerlab.github.io/syri/fileformat.html) for information on interpreting this table.
- `data/07_syri_output/<assembly>/<assembly>.plot.2.pdf` is a synteny plot of the rearrangements in the above tsv file. 

This command also generates other files. Please see the [output](output.md) page. 

>Tip: The 4th and 5th columns of the `<assembly>_clean.syri.out` file have the reference/query sequences for indels. For large deletions, these field can have thousands of bases. To view this tsv file in the terminal without these columns, try running this command. Note that older versions of the UNIX `cut` command may not support the `--complement` option, and you may need to either update `cut` or use other UNIX commands to hide the 4th and 5th columns. 

```
cat <assembly>_clean.syri.out | cut --complement -f 4-5
```

## Annotate genes in SVs

This command uses [prokka](https://github.com/tseemann/prokka) annotate the assemblies and [breseq](https://github.com/barricklab/breseq) to annotate the genes present in the deletions, inversions and amplifications.

```
snakemake --use-conda --cores n annotate_SV_regions
```

For each assembly, an HTML table `data/12_genome_diff_tables/html/<assembly>.html` . Please see the [_breseq_ manual]() for information on interpreting this table. For a sample table, please see [the tutorial](tutorial.md). 

## Predict SV mechanism

Most structural variant mutations occur through recombination between homologous sites and insertion sequences in particular are often involved. This command annotates the insertion sequences at the boundaries of structural variants (for deletions and inversions only) and predicts putative mechanisms.

```
snakamake --use-conda --cores n predict_SV_mechanism
```

- This command generates two summary csv files: `data/11_annotated_boundaries/deletion_mechanism.csv` and `data/11_annotated_boundaries/inversion_mechanism.csv` that contain information about the mechanism of deletions and inversions respectively, for all of the assemblies.
- Detailed information about the mechanism of each deletion for every assembly is in the csv files `data/11_annotated_boundaries/<assembly>_deletion.csv`
- Detailed information about the mechanism of each inversion for every assembly is in the csv files `data/11_annotated_boundaries/<assembly>_inversion.csv`

Information about these csv files can be found on the [output](output.md) page and examples on the [tutorial](tutorial.md) page.

## Predict replichore and inversion balance

For bacterial genomes that have a single origin and terminus, the genome can be divided into two replichores (halves) demarcated by the origin-terminus axis. It can be useful to know how the lengths of these two replichores differ between the ancestor-assembly pairs. The origin-terminus axis can also be used to classify inversions as inter-replichore or intra-replichore, and the symmetry of inter-replichore inversions about the origin-terminus axis can be measured. For an illustration explaining this, please see the [tutorial](tutorial.md).

This step requires user input (in addition to `data/data.csv`) to specify the sequences of the origin and terminus. For the origin, we recommend the *oriC* sequence, and for the terminus, we recommend the *dif* sequence. _seabreeze_ looks for these sequences (or the reverse complement) in the genome for an exact and unique match, and will cause an error if an exact match is not found, or more than one exact match is found. There is no minimum/maximum length requirements for the supplied _oriC_ and _dif_ sequences as long as they meet the above criteria. _seabreeze_ requires these two sequences for each ancestor, and assumes that the sequences have not mutated in the corresponding ancestor. This information is specified in a csv file in `data/ori_dif_sequences.csv`. This csv file has three columns: 

- `ancestor`: name of the ancestor, which matches the name in `data/data.csv`. It should not contain the `.fasta` extension.
- `ori`: base sequence of the _oriC_ locus, in upper case with only the standard bases (A,T,G,C). 
- `dif`: base sequence of the _dif_ locus, in upper case with only the standard bases (A,T,G,C). 

For example, if this is the `data/data.csv`:

| assembly | ancestor |
| -------- | -------- |
| genome1  | genome2  |
| genome3  | genome4  |

Then `data/ori_dif_sequence.csv` should look like this, with one entry per ancestor. The order of these rows does not matter, and each ancestor only requires one entry, even if it appears multiple times in `data/data.csv`. 

| ancestor | ori                | dif                |
| -------- | ------------------ | ------------------ |
| genome2  | GGATCCTGGGTATTAAAA | TCTTCCTTGGTTTATATT |
| genome4  | GTAGGATCCGTGATTAG  | AATCTGTCTCTGCACGTA |
With this csv file set up, this command can be run from the _seabreeze_ root directory:

```
snakemake --use-conda --cores n predict_replichore_balance
```

This command produces the following output:
- `data/08_reindex_genome_oric/replichore_arms.csv`: a csv file that describes the length of the two replichores for all of the assemblies
- `data/11_annotated_boundaries/inversion_replichores.csv`: a csv file that tallies the count of each type of inversion (inter-replichore vs intra-replichore) for all of the assemblies.
- `data/11_annotated_boundaries/<assembly>_inversion_classification.csv`: a csv file that has more detailed information about all of the inversions in a given assembly
- `data/11_annotated_boundaries/inversion_replichores_long.csv`: a csv file that contains a list of all inversions across all assemblies

Please see the [output](output.md) page for more information about these fields. 

---

## Additional Snakemake options

These are additional Snakemake options that can be included as follows. They work for all commands

```
snakemake --use-conda --cores n <insert flag here> <command>
```

- `-np`: This flag performs a 'dry-run' of the command i.e. it checks to see that all of the input files necessary for that command are present, or can be generated by _seabreeze_. It can helpful to run this before a computationally expensive step such as predicting IS elements, or predicting SVs.

If a command was terminated before it was completed, Snakemake may prompt you to use these flags when you run the same command again:

- `--unlock`:  Snakemake will "lock" a directory if all of the output files were not completed before the processes were terminated. Use this flag to unlock the directory.
- `--rerun-incomplete`:  To allow Snakemake to remake any output files that were not completed in the prior run.
