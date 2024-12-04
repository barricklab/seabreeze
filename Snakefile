import pandas as pd

path_to_data_csv="data/data.csv"
df = pd.read_csv(path_to_data_csv)

# TODO: run QC on the format of the data.csv

# this dictionary maps the subject to its query
assembly_to_ancestor_dict = dict(zip(df['assembly'], df['ancestor']))

# one rule to rule them all ...
# remember you cant have wildcards in the target rule!
rule all_targets:
    input:
        genome_sizes = "data/04_rename_genome/genome_size_stats.csv",
        #IS_summary = "data/05_isescan_tables/IS_summary.csv",
        #IS_summary_copy_change = "data/05_isescan_tables/IS_summary_copy_change.csv",
        #inversion_replichores = expand("data/11_annotated_boundaries/{sample}_inversion_classification.csv", sample=df['assembly'].tolist()),
        clean_synteny_plots = expand("data/07_syri_output/{sample}/{sample}.plot.2.pdf", sample=df['assembly'].tolist()),
        ori_dif_coords = "data/04_rename_genome/ori_dif_coords.csv",
        #ori_dif_coords = "data/08_reindex_genome_oric/ori_dif_coords.tsv",
        #replichore_arms = "data/08_reindex_genome_oric/replichore_arms.tsv",
        #deletion = expand("data/11_annotated_boundaries/{sample}_deletion.csv",sample=df['assembly'].tolist()),
        inversion_table = "data/11_annotated_boundaries/inversion_mechanism.csv",
        deletion_table = "data/11_annotated_boundaries/deletion_mechanism.csv",
        #inversion_classification = expand("data/11_annotated_boundaries/{sample}_inversion_classification.csv",sample=df['assembly'].tolist()),
        gd = expand("data/12_genome_diff_tables/gd/{sample}.gd",sample=df['assembly'].tolist()),
        html = expand("data/12_genome_diff_tables/html/{sample}.html",sample=df['assembly'].tolist())

# find unique bases at the start of the subject sequence to reindex the query sequence tp
rule find_reindex_bases:
    conda:
        "bin/workflow/envs/biopython.yml"
    input:
        query_path = "data/02_genomes/{sample}.fasta", # path to the assembly
        subject_path = lambda wildcards: "data/02_genomes/{}.fasta".format(assembly_to_ancestor_dict[wildcards.sample]), # path to the assembly its being compared to
        script = "bin/scripts/find_reindex_bases.py"
    output:
        "data/03_reindex_genomes/reindex_bases_{sample}.txt"
    log:
        "data/logs/find_reindex_bases/{sample}.log"
    shell:
        '''
        {input.script} --subject {input.subject_path} --query {input.query_path} --output {output} > {log} 2>&1
        '''

#reindex all the fasta file to a common sequence to make comparison easier
rule reindex_contigs:
    conda:
        "bin/workflow/envs/biopython.yml"
    input:
        fasta = "data/02_genomes/{sample}.fasta",
        bases = "data/03_reindex_genomes/reindex_bases_{sample}.txt",
        script = "bin/scripts/reindex_assembly.py"
    output:
       "data/03_reindex_genomes/{sample}.fasta"
    log:
        "data/logs/reindex_contigs/{sample}.log"
    shell:
       "{input.script} -b $(cat {input.bases})  -i {input.fasta} -o {output} -t fasta > {log} 2>&1"

# rename all the FASTA headers to "genome"
# this step is needed for SyRI which will only carry out variant calling for two sequence with the same header

rule rename_contigs:
    conda:
        "bin/workflow/envs/biopython.yml"
    input:
        data = "data/03_reindex_genomes/{sample}.fasta",
        script = "bin/scripts/rename_contigs.py"
    params:
        new_FASTA_header = "genome"
    output:
        "data/04_rename_genome/{sample}.fasta"
    log:
        "data/logs/rename_contigs/{sample}.log"
    shell:
        "{input.script} --file {input.data}  --name {params.new_FASTA_header} --output {output} > {log} 2>&1"


rule compute_genome_stats:
    conda:
        "bin/workflow/envs/biopython.yml"
    input:
        data =  expand("data/04_rename_genome/{sample}.fasta", sample=df['assembly'].tolist()),
        script = "bin/scripts/fasta_stats.py",
        csv_file = "data/data.csv"
    params:
        folder = "data/04_rename_genome"
    output:
        genome_sizes = "data/04_rename_genome/genome_size_stats.csv"
    log:
        "data/logs/compute_genome_stats/compute_genome_stats.log"
    shell:
        "{input.script} --folder {params.folder} --data {input.csv_file} --output {output.genome_sizes} > {log} 2>&1"

# ISEScan takes the genome assemblies and returns several files. We only need to the csv file it generates

rule find_IS_elements:
    conda:
        "bin/workflow/envs/isescan.yml"
    input:
        #expand("data/04_rename_genome/{sample}.fasta", sample=df['assembly'].tolist()) # dynamically generate list from csv file, not with wildcards
        "data/04_rename_genome/{sample}.fasta"
    output:
        #folder = "data/05_isescan_tables",
        csv_files = "data/05_isescan_tables/{sample}.csv"
        #csv_files = expand("data/05_isescan_tables/{sample}.csv", sample=df['assembly'].tolist()) # we only want the csv file, so that's the target of this rule is that.
    log:
        "data/logs/find_IS_elements/{sample}.log"
    shell:
        """
        echo {wildcards.sample}
        cp {input} ./{wildcards.sample}.fasta
        isescan.py --seqfile {wildcards.sample}.fasta --output data/05_isescan_tables/{wildcards.sample} --nthread 4 > {log} 2>&1
        mv data/05_isescan_tables/{wildcards.sample}/{wildcards.sample}.fasta.csv data/05_isescan_tables/{wildcards.sample}.csv
        rm {wildcards.sample}.fasta
        """


# align each assembly to its ancestor, then filter the alignments and convert from .delta to coords

rule align_genomes_nucmer:
    conda:
        "bin/workflow/envs/mummer4.yml"
    input:
        query_path = "data/04_rename_genome/{sample}.fasta", # path to the assembly
        subject_path = lambda wildcards: "data/04_rename_genome/{}.fasta".format(assembly_to_ancestor_dict[wildcards.sample]) # path to the assembly of the ancestor its being compared to
    output:
        done = "data/06_nucmer_alignment/{sample}/{sample}.done",
        delta = "data/06_nucmer_alignment/{sample}/{sample}.delta",
        filtered = "data/06_nucmer_alignment/{sample}/{sample}.filtered.delta",
        coords = "data/06_nucmer_alignment/{sample}/{sample}.filtered.coords"
    params:
        seq_id_cutoff = "95",
        subject_name = lambda wildcards: "{}.fasta".format(assembly_to_ancestor_dict[wildcards.sample]), # just the name of the ancestor (does not include the .fasta extension)
        output_dir = "data/06_nucmer_alignment/{sample}" # each assembly gets its own directory with the same name which stores the output of nucmer
    log:
        "data/logs/align_genomes_nucmer/{sample}.log"
    # temporarily move both fasta files here because it's easier. delete when done. NO DO NOT DO THIS! CAUSES A BUG WHEN COMPARING A SEQUENCE TO ITSELF
    shell:
        """
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        touch ../../../{log}
        nucmer --maxmatch -c 100 -b 500 -l 50 -p {wildcards.sample} ../../../{input.subject_path} ../../../{input.query_path} > ../../../{log} 2>&1
        delta-filter -i {params.seq_id_cutoff} -l 100 {wildcards.sample}.delta > {wildcards.sample}.filtered.delta
        show-coords -THrd {wildcards.sample}.filtered.delta > {wildcards.sample}.filtered.coords
        touch {wildcards.sample}.done
        cd ../../../
        echo "Alignment complete. Working dir set to:"
        pwd
        """

# now call structural variants from the alignments

rule call_variants_syri:
    conda:
        "bin/workflow/envs/syri.yml"
    input:
        filtered = "data/06_nucmer_alignment/{sample}/{sample}.filtered.delta",
        query_path = "data/04_rename_genome/{sample}.fasta", # path to the assembly
        subject_path = lambda wildcards: "data/04_rename_genome/{}.fasta".format(assembly_to_ancestor_dict[wildcards.sample]), # path to the assembly of the ancestor its being compared to
        coords = "data/06_nucmer_alignment/{sample}/{sample}.filtered.coords"
    output:
        done = "data/07_syri_output/{sample}/{sample}.done",
        syri = "data/07_syri_output/{sample}/{sample}syri.out"
    params:
        output_dir = "data/07_syri_output/{sample}", # each assembly gets its own directory with the same name which stores the output of nucmer
    log:
        "data/logs/call_variants_syri/{sample}.log"
    shell:
        """
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        touch ../../../{log}
        echo " the subject is {input.subject_path}"
        echo "the query is {input.query_path}"
        syri --nosnp -c ../../../{input.coords} -d ../../../{input.filtered} -r ../../../{input.subject_path} -q ../../../{input.query_path} --prefix {wildcards.sample} > ../../../{log} 2>&1
        touch {wildcards.sample}.done
        rm {wildcards.sample}syri.log {wildcards.sample}syri.summary
        cd ../../../
        echo "syri complete. Working dir set to:"
        pwd
        """

# generate the synteny plots with plotsr
# you start by creating the {sample}.genomes.tsv file needed by plotsr. This is created for each sample
rule generate_synteny_plot:
    conda:
        "bin/workflow/envs/plotsr.yml"
    input:
        query_path = "data/04_rename_genome/{sample}.fasta", # path to the assembly
        subject_path = lambda wildcards: "data/04_rename_genome/{}.fasta".format(assembly_to_ancestor_dict[wildcards.sample]), # path to the assembly of the ancestor its being compared to
        syri = "data/07_syri_output/{sample}/{sample}syri.out",
        script = "bin/scripts/plotsr/plotsr-bin"
    output:
        genome_table = "data/07_syri_output/{sample}/{sample}.genomes.tsv",
        plot = "data/07_syri_output/{sample}/{sample}.plot.pdf"
    params:
        input_dir = "data/07_syri_output/{sample}", #store the synteny plot in the same place as the syri files
        subject_name = lambda wildcards: "{}.fasta".format(assembly_to_ancestor_dict[wildcards.sample]) # just the name of the ancestor (does not include the .fasta extension)
    log:
        "data/logs/generate_synteny_plots/{sample}.log"
    shell:
        """
        cd {params.input_dir}
        printf "#file\tname\ttags\n" > {wildcards.sample}.genomes.tsv
        printf "../../../{input.subject_path}\t{params.subject_name}\tlw:1.5\n" >> {wildcards.sample}.genomes.tsv
        printf "../../../{input.query_path}\t{wildcards.sample}\tlw:1.5" >>  {wildcards.sample}.genomes.tsv
        ../../../{input.script} -s 500 --genomes {wildcards.sample}.genomes.tsv --sr {wildcards.sample}syri.out -H 5 -W 10 -o {wildcards.sample}.plot.pdf --lf {wildcards.sample}.log
        mv {wildcards.sample}.log ../../../{log}
        cd ../../..
        echo "Synteny plot generated. Working dir set to:"
        pwd
        """

# now clean up the syri files to predict a minimal set of structural variants
rule clean_syri_output:
    conda:
        "bin/workflow/envs/pandas.yml"
    input:
        syri = "data/07_syri_output/{sample}/{sample}syri.out",
        query_path = "data/05_isescan_tables/{sample}.csv", # path to the isescan file of the
        subject_path = lambda wildcards: "data/05_isescan_tables/{}.csv".format(assembly_to_ancestor_dict[wildcards.sample]) # path to the assembly of the ancestor its being compared to
    output:
        "data/07_syri_output/{sample}/{sample}syri.out_v2"
    params:
        # isescan_subject_path = expand("data/05_isescan_tables/{sample}.csv", sample=df['assembly'].tolist()), # listing this as an input triggers an InputExceptionError idk why
        # isescan_query = lambda wildcards: "{}.csv".format(assembly_to_ancestor_dict[wildcards.sample]), # just the name of the ancestor (does not include the .fasta extension)
        isescan_dir = "data/05_isescan_tables",
        input_dir = "data/07_syri_output/{sample}",
        script = "bin/scripts/clean_syri.py"
    log:
        "data/logs/clean_syri_output/{sample}.log"
    shell:
        """
        cd {params.input_dir}
        echo "{input.subject_path}"
        ../../../{params.script} --syri {wildcards.sample}syri.out --isescan_query ../../../{input.query_path} --isescan_subject ../../../{input.subject_path} > ../../../{log} 2>&1
        cd ../../..
        echo "working dir set back to"
        pwd
        """
# with the new clean syri file, generate a new plot
rule generate_synteny_plot_clean:
    conda:
        "bin/workflow/envs/plotsr.yml"
    input:
        syri = "data/07_syri_output/{sample}/{sample}syri.out_v2",
        script = "bin/scripts/plotsr/plotsr-bin",
        genome_table = "data/07_syri_output/{sample}/{sample}.genomes.tsv",
    output:
        "data/07_syri_output/{sample}/{sample}.plot.2.pdf"
    params:
        input_dir = "data/07_syri_output/{sample}", #store the synteny plot in the same place as the syri files
        subject_name = lambda wildcards: "{}.fasta".format(assembly_to_ancestor_dict[wildcards.sample]) # just the name of the ancestor (does not include the .fasta extension)
    log:
        "data/logs/generate_synteny_plots/{sample}.2.log"
    shell:
        """
        cd {params.input_dir}
        ../../../{input.script} -s 500 --genomes {wildcards.sample}.genomes.tsv --sr {wildcards.sample}syri.out_v2 -H 5 -W 10 -o {wildcards.sample}.plot.2.pdf --lf {wildcards.sample}.2.log
        mv {wildcards.sample}.2.log ../../../{log}
        cd ../../..
        pwd
        """

# generate a csv file with the ori and dif coords of the genomes in their original index
rule annotate_ori_dif_locations:
    conda:
        "bin/workflow/envs/pandas.yml"
    input:
        genomes = expand("data/04_rename_genome/{sample}.fasta", sample=df['assembly'].tolist()), # you can't use wildcards here but you can use this expand functionality
        script = "bin/scripts/replichore_arms_analyse.py"
    output:
        ori_dif_coords = "data/04_rename_genome/ori_dif_coords.csv"
    params:
        data="data/data.csv",
        sequences="data/sequences.csv",
        folder = "data/04_rename_genome/"
    shell:
        """
        {input.script} --genomes {params.folder} --data {params.data} --sequences {params.sequences} --output {output} --noarms
        pwd
        """

# reindex all the fasta file to the origin to analyse the replichore arms and find ori and dif position
# TODO: The input for this should be ori_dif_coords.tsv as each genome's ori location should be fetched from there

rule reindex_contigs_oric:
    conda:
        "bin/workflow/envs/biopython.yml"
    input:
        data = "data/02_genomes/{sample}.fasta",
        script = "bin/scripts/reindex_assembly.py"
    output:
        "data/08_reindex_genome_oric/{sample}.fasta"
    shell:
        "{input.script} -b GGATCCTGGGTATTAAAA -i {input.data} -o {output} -t fasta"


# generate a tsv file with the oric and dif of the genomens reindexed to the ori and a tsv file with lenths of the replichore arms of each clone
# TODO: Eventually remove the ancestor as a requirement for this rule and script, since it is not used

rule analyse_replichore_arms:
    conda:
        "bin/workflow/envs/pandas.yml"
    input:
        genomes = expand("data/08_reindex_genome_oric/{sample}.fasta", sample=df['assembly'].tolist()), # you can't use wildcards here but you can use this expand functionality
        script = "bin/scripts/replichore_arms_analyse.py"
    output:
        ori_dif_coords = "data/08_reindex_genome_oric/ori_dif_coords.tsv",
        replichore_arms = "data/08_reindex_genome_oric/replichore_arms.tsv"
    params:
        folder = "data/08_reindex_genome_oric/",
        ori = "GGATCCTGGGTATTAAAA",
        dif = "TCTTCCTTGGTTTATATT",
        ancestor = "Anc-_0gen_REL606",
        output_name_oridif = "ori_dif_coords.tsv",
        output_name_arms = "replichore_arms.tsv"

    shell:
        """
        {input.script} --assemblies {params.folder} --ori  {params.ori} --dif {params.dif} --ancestor {params.ancestor} --output {params.output_name_oridif} --noarms
        pwd
        {input.script} --assemblies {params.folder} --ori  {params.ori} --dif {params.dif} --ancestor {params.ancestor} --output {params.output_name_arms}
        """

# Run breseq to predict deletions and amplifications to see if they were missed
# i am deleting some of the output of breseq but feel free to remove that line in case you want it
# rule run_breseq:
#     conda:
#         "bin/workflow/envs/breseq.yml"
#     input:
#         reads = "data/09_merged_trimmed_nanopore_reads/{sample}.nanopore.fastq.gz", # reads of the evolved clone
#         reference_assembly = lambda wildcards: "data/04_rename_genome/{}.fasta".format(assembly_to_ancestor_dict[wildcards.sample])
#     output:
#         gd = "data/10_breseq_output/{sample}.gd",
#         html = "data/10_breseq_output/{sample}.html"
#     log:
#         "data/logs/run_breseq/{sample}.log"
#     params:
#         breseq_dir = "{sample}",
#         threads = "4",
#         limit_reads = "60" # this speeds up breseq by limiting the read depth to which it looks at data
#     shell:
#         """
#         mkdir -p data/10_breseq_output
#         cd data/10_breseq_output
#         breseq -j {params.threads} -l {params.limit_reads} -x -r ../../{input.reference_assembly} ../../{input.reads} -o {params.breseq_dir}> {wildcards.sample}.log 2>&1
#         mv {wildcards.sample}.log ../../{log}
#         cd {params.breseq_dir}
#         rm -rf 01_sequence_conversion 03_candidate_junctions 05_alignment_correction 07_error_calibration 02_reference_alignment 04_candidate_junction_alignment 06_bam 08_mutation_identification
#         mv data/output.gd ../{wildcards.sample}.gd
#         mv output/index.html ../{wildcards.sample}.html
#         cd ../../..
#         echo "task done. wd set to"
#         """

# generate tsv files which annotate the boundaries of the SVs
rule annotate_SV_boundaries_IS:
    conda:
        "bin/workflow/envs/pandas.yml"
    input:
        syri = "data/07_syri_output/{sample}/{sample}syri.out_v2",
        assembly_IS_csv = "data/05_isescan_tables/{sample}.csv",
        ancestor_IS_csv = lambda wildcards: "data/05_isescan_tables/{}.csv".format(assembly_to_ancestor_dict[wildcards.sample]), # path to the csv file of the ancestor
        script = "bin/scripts/IS_SV_border.py"
    output:
        "data/11_annotated_boundaries/{sample}_boundaries.tsv"
    shell:
        """
        mkdir -p data/11_annotated_boundaries
        cd data/11_annotated_boundaries
        ../../{input.script} --ancestor ../../{input.ancestor_IS_csv} --evolved ../../{input.assembly_IS_csv} --syri ../../{input.syri} --output {wildcards.sample}_boundaries.tsv
        cd ../..
        """

# annotate the mechanism of deletions and inversions
# i think you can expand and not wildcards here since the script does not have to be repeated each you run this..

rule annotate_SV_mechanism:
    conda:
        "bin/workflow/envs/biopython.yml"
    input:
        boundaries_csv = expand("data/11_annotated_boundaries/{sample}_boundaries.tsv", sample=df['assembly'].tolist()),
        script = "bin/scripts/classify_deletions.py"
    output:
        inversion = expand("data/11_annotated_boundaries/{sample}_inversion.csv",sample=df['assembly'].tolist()),
        deletion = expand("data/11_annotated_boundaries/{sample}_deletion.csv",sample=df['assembly'].tolist()),
        inversion_table = "data/11_annotated_boundaries/inversion_mechanism.csv",
        deletion_table = "data/11_annotated_boundaries/deletion_mechanism.csv"
    params:
        input_dir = "data/11_annotated_boundaries/",
        output_deletion = "deletion_mechanism.csv",
        output_inversion = "inversion_mechanism.csv"
    shell:
        """
        {input.script} --folder {params.input_dir} --output {params.output_inversion} --inversion
        {input.script} --folder {params.input_dir} --output {params.output_deletion} --deletion
        cd ..
        """

# classify inversions as inter_replichore or intra-replichore

rule classify_inversion_replichore:
    conda:
        "bin/workflow/envs/biopython.yml"
    input:
        ori_dif_coords = "data/04_rename_genome/ori_dif_coords.tsv",
        inversion = expand("data/11_annotated_boundaries/{sample}_inversion.csv",sample=df['assembly'].tolist()),
        #inversion = "data/11_annotated_boundaries/{sample}_inversion.csv",
        script = "bin/scripts/inversion_replichore_classify.py"
    output:
        #"data/11_annotated_boundaries/{sample}_inversion_classification.csv"
        expand("data/11_annotated_boundaries/{sample}_inversion_classification.csv",sample=df['assembly'].tolist())
    params:
        input_dir = "data/11_annotated_boundaries/",
        output_table = "inversion_replichores.csv",
    shell:
        """
        {input.script} --folder {params.input_dir} --oridif {input.ori_dif_coords} --output {params.output_table} --data data.csv
        """

# generate annotations for genomes using prokka and with the IS elements reported by ISEscan

rule annotate_genomes_prokka:
    conda:
        "bin/workflow/envs/prokka.yml"
    input:
        genome="data/04_rename_genome/{sample}.fasta",
        is_table="data/05_isescan_tables/{sample}.csv"
    output:
        "09_annotated_genomes/{sample}.gff"
    params:
        prefix="{sample}",
        outdir="09_annotated_genomes/{sample}",
        prokka_annotation="09_annotated_genomes/{sample}/{sample}.gff"
    log:
        "data/logs/annotate_genomes_prokka/{sample}.log"
    shell:
        """
        prokka --prefix {params.prefix} --outdir {params.outdir} {input.genome} > {log} 2>&1
        breseq CONVERT-REFERENCE -f GFF3 -s {input.is_table} -o {output} {params.prokka_annotation} >> {log} 2>&1
        """


# Use the syri.out_v2 files to make the HTML tables from breseq
rule generate_genome_diffs_tables:
    conda:
        "bin/workflow/envs/breseq.yml"
    input:
        syri = "data/07_syri_output/{sample}/{sample}syri.out_v2",
        script = "bin/scripts/syri2gd.py",
        reference = "09_annotated_genomes/{sample}.gff"
    output:
        gd = "data/12_genome_diff_tables/gd/{sample}.gd",
        html = "data/12_genome_diff_tables/html/{sample}.html"
    params:
        gd_folder = "data/12_genome_diff_tables/gd",
        html_folder = "data/12_genome_diff_tables/html",
    shell:
        """
        mkdir -p {params.gd_folder}
        mkdir -p {params.html_folder}
        cd {params.gd_folder}
        ../../../{input.script} --syri ../../../{input.syri} --output {wildcards.sample}.gd --deletion --inversion --amplification
        cd ../../../{params.html_folder}
        gdtools ANNOTATE -o {wildcards.sample}.html -r ../../../{input.reference} -f HTML ../../../{output.gd}
        cd ../../../
        """
