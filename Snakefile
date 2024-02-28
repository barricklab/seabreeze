import pandas as pd

path_to_csv="data.csv"
df = pd.read_csv(path_to_csv)

# these dictionaries map assembly to its path and assembly to its ancestor
assembly_to_path_dict = dict(zip(df['assembly'], df['assembly_path']))
assembly_to_ancestor_dict = dict(zip(df['assembly'], df['ancestor']))


# reindex all the fasta file to a common sequence to make comparison easier
rule reindex_contigs:
    conda:
        "bin/workflow/envs/biopython.yml"
    input:
        data = "data/02_genomes/{sample}.fasta",
        script = "bin/scripts/reindex_assembly.py"
    output:
        "data/03_reindex_genome/{sample}.fasta"
    shell:
        "{input.script} -b AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTC -i {input.data} -o {output} -t fasta"

# rename all the contigs of the fasta files to a string (here "REL606")
# this step is needed for SyRI which will only carry out variant calling for two sequence with the same header
rule rename_contigs:
    conda:
        "bin/workflow/envs/biopython.yml"
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
        "bin/workflow/envs/pandas.yml"
    input:
        data = "data/04_rename_genome",
        clones = expand("data/05_isescan_tables/{sample}.csv", sample=df['assembly'].tolist()),
        script = "bin/scripts/fasta_stats.py"
    params:
        ancestor = "Anc-_0gen_REL606.fasta" # this is not the path to the ancestor's assembly but it is expected that the ancestor is in input.data folder
    output:
        contig_stats = "data/04_rename_genome/contig_stats.tsv",
        genome_sizes = "data/04_rename_genome/genome_size_stats.tsv"
    shell:
        "{input.script} --folder {input.data} --output {output.contig_stats} --stats {output.genome_sizes} --ancestor {params.ancestor}"
 
 #ISEScan takes the genome assemblies and returns several files. We only need to the csv file it generates
rule find_IS_elements:
    conda:
        "bin/workflow/envs/isescan.yml"
    input:
        "data/04_rename_genome/{sample}.fasta"
    output:
        "data/05_isescan_tables/{sample}.csv" # we only want the csv file, so that's the target of this rule is that. 
    shell:
        """
        echo {wildcards.sample}
        cp {input} ./{wildcards.sample}.fasta
        isescan.py --seqfile {wildcards.sample}.fasta --output data/05_isescan_tables/{wildcards.sample} --nthread 4
        mv data/05_isescan_tables/{wildcards.sample}/{wildcards.sample}.fasta.csv data/05_isescan_tables/{wildcards.sample}.csv
        rm {wildcards.sample}.fasta
        """

# align each assembly to its ancestor, then filter the alignments and convert from .delta to coords
rule align_genomes_nucmer:
    conda:
        "bin/workflow/envs/mummer4.yml"
    input:
        subject_path = "data/04_rename_genome/{sample}.fasta", # path to the assembly 
        query_path = lambda wildcards: "data/04_rename_genome/{}.fasta".format(assembly_to_ancestor_dict[wildcards.sample]) # path to the assembly of the ancestor its being compared to
    output:
        done = "data/06_nucmer_alignment/{sample}/{sample}.done",
        delta = "data/06_nucmer_alignment/{sample}/{sample}.delta",
        filtered = "data/06_nucmer_alignment/{sample}/{sample}.filtered.delta",
        coords = "data/06_nucmer_alignment/{sample}/{sample}.filtered.coords"
    params:
        seq_id_cutoff = "95",
        query_name = lambda wildcards: "{}.fasta".format(assembly_to_ancestor_dict[wildcards.sample]), # just the name of the ancestor (does not include the .fasta extension)
        output_dir = "data/06_nucmer_alignment/{sample}", # each assembly gets its own directory with the same name which stores the output of nucmer
    log:
        "data/logs/align_genomes_nucmer/{sample}.log"
# temporarily move both fasta files here because it's easier. delete when done.
    shell:
        """
        cp {input.subject_path} {params.output_dir}/{wildcards.sample}.fasta 
        cp {input.query_path} {params.output_dir}/{params.query_name}.fasta
        cd data/06_nucmer_alignment/{wildcards.sample}
        nucmer --maxmatch -c 100 -b 500 -l 50 -p {wildcards.sample} {params.query_name}.fasta {wildcards.sample}.fasta > {log}
        delta-filter -i {params.seq_id_cutoff} -l 100 {wildcards.sample}.delta > {wildcards.sample}.filtered.delta
        show-coords -THrd {wildcards.sample}.filtered.delta > {wildcards.sample}.filtered.coords
        touch {wildcards.sample}.done
        rm {wildcards.sample}.fasta {params.query_name}.fasta
        cd ../../..
        """

# TODO: The target rule generates o/p files based on the data.csv file. Use a function which returns a list
# rule all:
#     input:
#         expand("data/04_rename_genome/{output_tsv}", output_tsv=["contig_stats.tsv", "genome_size_stats.tsv"]),
#         expand("data/05_isescan_tables/{sample}.csv", sample=df['clone'].tolist())
