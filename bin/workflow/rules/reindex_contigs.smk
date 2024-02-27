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
