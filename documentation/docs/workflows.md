General usage for a _seabreeze_ command can be found on [usage](usage.md). All output files and folder and generated in the specified working directory. Depending on if the `--mask` option was enabled, output will be stored in subdirectories called `masked` or `unmasked`. For clarity, this part of the output file path has been omitted in this documentation, but the users should find navigating the output intuitive nonetheless.

## Run all

Runs all of the following workflows in a single command.

## Analyse genome sizes

This workflow compares the sizes of the assemblies to their specified ancestors in a pairwise manner. 

The output file generated is `04_rename_genome/genome_size_stats.csv`. Please see [output](output.md) for more information about the fields in this table, and the [tutorial](tutorial.md) for a sample.
## Find insertion sequences

This workflow uses [ISEscan](https://github.com/xiezhq/ISEScan/blob/master/README.md?plain=1) to annotate partial and full-length IS elements in the genomes. It primarily generates a series of csv files`05_isescan_tables/<assembly>.csv`, with each csv file corresponding to an assembly in `data.csv` (if running in batch mode). ISEScan also generates additional output for each assembly, which is stored in their respective directories `05_iescan_tables/<assembly>/`.  For information about interpreting this csv file, we refer users to the official [ISEScan documentation](https://github.com/xiezhq/ISEScan/blob/master/README.md). 
## Predict SVs

This command aligns the assembly-ancestor pairs, predicts structural variants and generates synteny plots. 
- `07_syri_output/<assembly>/<assembly>_clean.syri.out` is a tsv file that describes the location and type of SVs in the assembly, relative to its ancestor. We refer the user to the official [SyRI documentation](https://schneebergerlab.github.io/syri/fileformat.html) for information on interpreting this table.
- `07_syri_output/<assembly>/<assembly>_clean.plot.pdf` is a synteny plot of the rearrangements in the above tsv file. 

This command also generates other files. Please see the [output](output.md) page. 

>Tip: The 4th and 5th columns of the `<assembly>_clean.syri.out` file have the reference/query sequences for indels. For large deletions, these field can have thousands of bases. To view this tsv file in the terminal without these columns, try running this command. Note that older versions of the UNIX `cut` command may not support the `--complement` option, and you may need to either update `cut` or use other UNIX commands to hide the 4th and 5th columns. 

```
cat <assembly>_clean.syri.out | cut --complement -f 4-5
```

## Annotate genes in SVs

This command uses [prokka](https://github.com/tseemann/prokka) annotate the assemblies and [breseq](https://github.com/barricklab/breseq) to annotate the genes present in the deletions, inversions and amplifications.

For each assembly, an HTML table `12_genome_diff_tables/html/<assembly>.html` . Please see the [_breseq_ manual]() for information on interpreting this table. For a sample table, please see [the tutorial](tutorial.md). 

## Predict SV mechanism

Most structural variant mutations occur through recombination between homologous sites and insertion sequences in particular are often involved. This command annotates the insertion sequences at the boundaries of structural variants (for deletions and inversions only) and predicts putative mechanisms.

- This command generates two summary csv files: `11_annotated_boundaries/deletion_mechanism.csv` and `11_annotated_boundaries/inversion_mechanism.csv` that contain information about the mechanism of deletions and inversions respectively, for all of the assemblies.
- Detailed information about the mechanism of each deletion for every assembly is in the csv files `11_annotated_boundaries/<assembly>_deletion.csv`
- Detailed information about the mechanism of each inversion for every assembly is in the csv files `11_annotated_boundaries/<assembly>_inversion.csv`

Information about these csv files can be found on the [output](output.md) page and examples on the [tutorial](tutorial.md) page.

## Predict replichore and inversion balance

For bacterial genomes that have a single origin and terminus, the genome can be divided into two replichores (halves) demarcated by the origin-terminus axis. It can be useful to know how the lengths of these two replichores differ between the ancestor-assembly pairs. The origin-terminus axis can also be used to classify inversions as inter-replichore or intra-replichore, and the symmetry of inter-replichore inversions about the origin-terminus axis can be measured. For an illustration explaining this, please see the [tutorial](tutorial.md).

This step requires user input (in addition to `data.csv`) to specify the sequences of the origin and terminus. For the origin, we recommend the *oriC* sequence, and for the terminus, we recommend the *dif* sequence. _seabreeze_ looks for these sequences (or the reverse complement) in the genome for an exact and unique match, and will cause an error if an exact match is not found, or more than one exact match is found. There is no minimum/maximum length requirements for the supplied _oriC_ and _dif_ sequences as long as they meet the above criteria. It also assumes that the sequences have not mutated in the corresponding assembly.

If using the `seabreeze run` command, these are supplied with the `--ori` and `--dif` options.

If using the `seabreeze batch` command, it needs these two sequences for each ancestor. This information is specified in a csv file in `data/ori_dif_sequences.csv`. This csv file has three columns: 

- `ancestor`: name of the ancestor, which matches the name in `data.csv`. It should not contain the `.fasta` extension.
- `ori`: base sequence of the _oriC_ locus,
- `dif`: base sequence of the _dif_ locus

For example, if this is the `data.csv`:

| assembly | ancestor |
| -------- | -------- |
| genome1  | genome2  |
| genome3  | genome4  |

Then `ori_dif_sequence.csv` should look like this, with one entry per ancestor. The order of these rows does not matter, and each ancestor only requires one entry, even if it appears multiple times in `data.csv`. 

| ancestor | ori                | dif                |
| -------- | ------------------ | ------------------ |
| genome2  | GGATCCTGGGTATTAAAA | TCTTCCTTGGTTTATATT |
| genome4  | GTAGGATCCGTGATTAG  | AATCTGTCTCTGCACGTA |

This command produces the following output:
- `08_reindex_genome_oric/replichore_arms.csv`: a csv file that describes the length of the two replichores for all of the assemblies
- `11_annotated_boundaries/inversion_replichores.csv`: a csv file that tallies the count of each type of inversion (inter-replichore vs intra-replichore) for all of the assemblies.
- `11_annotated_boundaries/<assembly>_inversion_classification.csv`: a csv file that has more detailed information about all of the inversions in a given assembly
- `11_annotated_boundaries/inversion_replichores_long.csv`: a csv file that contains a list of all inversions across all assemblies

Please see the [output](output.md) page for more information about these fields. 

---
