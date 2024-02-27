#ISEScan takes the genome assemblies and returns several files. We only need to the csv file it generates
rule find_IS_elements:
    conda:
        "bin/workflow/envs/isescan.yml"
    input:
        data = "data/04_rename_genome/{sample}.fasta"
    output:
        "data/05_isescan_tables/{sample}.csv"
    shell:
        """
        isescan.py --seqfile {input.data} --output {wildcards.sample} --nthread 4
        mv data/05_isescan_tables/{wildcards.sample}/{wildcards.sample}.fasta.csv data/05_isescan_tables/{wildcards.sample}.csv
        """