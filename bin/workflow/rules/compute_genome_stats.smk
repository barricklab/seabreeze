
# Calculate the number of contigs in each fasta file and their length. Output is stored in contig_stats.tsv
# Calculate the difference in length of the genomes, relative to the ancestor genome_size_stats.tsv

rule compute_genome_stats:
    conda:
        "bin/workflow/envs/pandas.yml"
    input:
        data = "data/04_rename_genome",
        clones = "data/04_rename_genome/{sample}.fasta"
        script = "bin/scripts/fasta_stats.py"
    params:
        ancestor = "Anc-_0gen_REL606.fasta" # this is not the path to the ancestor's assembly but it is expected that the ancestor is in input.data folder
    output:
        contig_stats = "data/04_rename_genome/contig_stats.tsv",
        genome_sizes = "data/04_rename_genome/genome_size_stats.tsv"
    shell:
        "{input.script} --folder {input.data} --output {output.contig_stats} --stats {output.genome_sizes} --ancestor {params.ancestor}"
