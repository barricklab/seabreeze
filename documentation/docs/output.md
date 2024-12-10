# Output

## Analyse genome sizes

This command generates `data/04_rename_genome/genome_size_stats.csv` as output. This csv file has the following fields:
assembly: name of the assembly, as specified in `data/data.csv`
ancestor: name of the ancestor, as specified in `data/data.csv`
size_assembly: size of the assembly, in bp
size_ancestor: size of the ancestor, in bp
difference: size_assembly - size_ancestor. This value is negative
percent_change: the difference in size, as a percent of the size of the assembly 

Here is a sample table generate from the [tutorial](tutorial.md)

|assembly|ancestor|size_assembly|size_ancestor|difference|percent_change|
|-----|-------|------|-------|----|----|
|REL606|REL606|4629812|4629812|0|0.0|
|REL606_evolved_1|REL606|4617111|4629812|-12701|-0.2743|
|REL606_evolved_2|REL606|4549910|4629812|-79902|-1.7258|
|REL606_evolved_3|REL606|4673209|4629812|43397|0.9373|

