# reindex all the fasta file to a common sequence to make comparison easier
rule reindex_contigs:
    conda:
        "bin/envs/biopython.yml"
    input:
        data = "data/02_genomes/{sample}.fasta",
        script = "bin/scripts/reindex_assembly.py"
    output:
        "data/03_reindex_genome/{sample}.fasta"
    shell:
        "{input.script} -b AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTC -i {input.data} -o {output} -t fasta"

# TODO: wrap this up to accept the name to reindex to
# rename all the contigs of the fasta files to a string (here "REL606")
# this step is needed for SyRI which will only carry out variant calling for two sequence with the same header
rule rename_contigs:
    conda:
        "bin/envs/biopython.yml"
    input:
        data = "data/03_reindex_genome/{sample}.fasta",
        script = "bin/scripts/rename_contigs.py"
    output:
        "data/04_rename_genome/{sample}.fasta"
    shell:
        "{input.script} --file --fasta {input.data}  --name REL606 --output {output}"

# Calculate the number of contigs in each fasta file and their length. Output is stored in contig_stats.tsv
# Calculate the difference in length of the genomes, relative to the ancestor genome_size_stats.tsv

rule compute_genome_stats:
    conda:
        "bin/envs/pandas.yml"
    input:
        data = "data/04_rename_genome",
        script = "bin/scripts/fasta_stats.py"
    params:
        ancestor = "Anc-_0gen_REL606.fasta" # this is not the path to the ancestor's assembly but it is expected that the ancestor is in input.data folder
    output:
        contig_stats = "data/04_rename_genome/contig_stats.tsv",
        genome_sizes = "data/04_rename_genome/genome_size_stats.tsv"
    shell:
        "{input.script} --folder {input.data} --output {output.contig_stats} --stats {output.genome_sizes} --ancestor {params.ancestor}"
