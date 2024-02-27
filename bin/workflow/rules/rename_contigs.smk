# TODO: wrap this up to accept the name to reindex to
# rename all the contigs of the fasta files to a string (here "REL606")
# this step is needed for SyRI which will only carry out variant calling for two sequence with the same header
rule rename_contigs:
    conda:
        "bin//workflow/envs/biopython.yml"
    input:
        data = "data/03_reindex_genome/{sample}.fasta",
        script = "bin/scripts/rename_contigs.py"
    output:
        "data/04_rename_genome/{sample}.fasta"
    shell:
        "{input.script} --file --fasta {input.data}  --name REL606 --output {output}"
