import pandas as pd

path_to_csv="data.csv"
df = pd.read_csv(path_to_csv)

# these dictionaries map assembly to its path and assembly to its ancestor
assembly_to_path_dict = dict(zip(df['assembly'], df['assembly_path']))
assembly_to_ancestor_dict = dict(zip(df['assembly'], df['ancestor']))

# one rule to rule them all ...
# remember you cant have wildcards in the target rule!
rule all_targets:
    input:
        contig_stats = "data/04_rename_genome/contig_stats.tsv",
        genome_sizes = "data/04_rename_genome/genome_size_stats.tsv",
        IS_summary = "data/05_isescan_tables/IS_summary.csv",
        IS_summary_copy_change = "data/05_isescan_tables/IS_summary_copy_change.csv",
        inversion_replichores = expand("data/11_annotated_boundaries/{sample}_inversion_classification.csv", sample=df['assembly'].tolist()),
        clean_synteny_plots = expand("data/07_syri_output/{sample}/{sample}.plot.2.pdf", sample=df['assembly'].tolist()),
        ori_dif_coords = "data/08_reindex_genome_oric/ori_dif_coords.tsv",
        replichore_arms = "data/08_reindex_genome_oric/replichore_arms.tsv",
        deletion = expand("data/11_annotated_boundaries/{sample}_deletion.csv",sample=df['assembly'].tolist()),
        inversion_table = "data/11_annotated_boundaries/inversion_mechanism.csv",
        deletion_table = "data/11_annotated_boundaries/deletion_mechanism.csv",
        inversion_classification = expand("data/11_annotated_boundaries/{sample}_inversion_classification.csv",sample=df['assembly'].tolist())

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
        data =  expand("data/04_rename_genome/{sample}.fasta", sample=df['assembly'].tolist()),
        script = "bin/scripts/fasta_stats.py"
    params:
        folder = "data/04_rename_genome",
        ancestor = "Anc-_0gen_REL606.fasta" # this is not the path to the ancestor's assembly but it is expected that the ancestor is in input.data folder
    output:
        contig_stats = "data/04_rename_genome/contig_stats.tsv",
        genome_sizes = "data/04_rename_genome/genome_size_stats.tsv"
    shell:
        "{input.script} --folder {params.folder} --output {output.contig_stats} --stats {output.genome_sizes} --ancestor {params.ancestor}"
 
 #ISEScan takes the genome assemblies and returns several files. We only need to the csv file it generates
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

# from the ISEScan tables, generate a summary of the total copy number, and the change in copy number relative to oone ancestor
# TODO: Update this so each assembly's copy number change is calculate for it's specified ancestor
rule generate_ISEScan_summary:
    conda:
        "bin/workflow/envs/pandas.yml"
    input:
        is_csv = expand("data/05_isescan_tables/{sample}.csv", sample=df['assembly'].tolist()), # you can't use wildcards here but you can use this expand functionality
        script = "bin/scripts/isescan_summary.py"
    output:
        "data/05_isescan_tables/IS_summary.csv",
        "data/05_isescan_tables/IS_summary_copy_change.csv"
    params:
        output_name = "IS_summary",
        ancestor = "Anc-_0gen_REL606",
        input_folder = "data/05_isescan_tables",
    shell:
        """
        {input.script} --isescan {params.input_folder} --output {params.output_name} --ancestor {params.ancestor}
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
        echo "thhe query is {input.query_path}"
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

# reindex all the fasta file to the origin to analyse the replichore arms and find ori and dif position
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

# generate a tsv file with the ori and dif coords 
rule annotate_ori_dif_locations:
    conda:
        "bin/workflow/envs/pandas.yml"
    input:
        genomes = expand("data/04_rename_genome/{sample}.fasta", sample=df['assembly'].tolist()), # you can't use wildcards here but you can use this expand functionality
        script = "bin/scripts/replichore_arms_analyse.py"
    output:
        ori_dif_coords = "data/04_rename_genome/ori_dif_coords.tsv"
    params:
        folder = "data/04_rename_genome/",
        ori = "GGATCCTGGGTATTAAAA",
        dif = "TCTTCCTTGGTTTATATT",
        ancestor = "Anc-_0gen_REL606",
        output_name_oridif = "ori_dif_coords.tsv"
    shell:
        """
        {input.script} --assemblies {params.folder} --ori  {params.ori} --dif {params.dif} --ancestor {params.ancestor} --output {params.output_name_oridif} --noarms
        pwd
        """

# generate a tsv file with the oric and dif and a tsv file with lenths of the replichore arms of each clone
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
rule run_breseq:
    conda:
        "bin/workflow/envs/breseq.yml"
    input:
        reads = "data/09_merged_trimmed_nanopore_reads/{sample}.nanopore.fastq.gz", # reads of the evolved clone
        reference_assembly = lambda wildcards: "data/04_rename_genome/{}.fasta".format(assembly_to_ancestor_dict[wildcards.sample])
    output:
        gd = "data/10_breseq_output/{sample}.gd",
        html = "data/10_breseq_output/{sample}.html"
    log:
        "data/logs/run_breseq/{sample}.log"
    params:
        breseq_dir = "{sample}",
        threads = "4",
        limit_reads = "60" # this speeds up breseq by limiting the read depth to which it looks at data
    shell:
        """
        mkdir -p data/10_breseq_output
        cd data/10_breseq_output
        breseq -j {params.threads} -l {params.limit_reads} -x -r ../../{input.reference_assembly} ../../{input.reads} -o {params.breseq_dir}> {wildcards.sample}.log 2>&1
        mv {wildcards.sample}.log ../../{log}
        cd {params.breseq_dir}
        rm -rf 01_sequence_conversion 03_candidate_junctions 05_alignment_correction 07_error_calibration 02_reference_alignment 04_candidate_junction_alignment 06_bam 08_mutation_identification 
        mv data/output.gd ../{wildcards.sample}.gd
        mv output/index.html ../{wildcards.sample}.html
        cd ../../..
        echo "task done. wd set to"
        """

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
rule annotate_SV_mechanism:
    conda:
        "bin/workflow/envs/biopython.yml"
    input:
        boundaries_csv = "data/11_annotated_boundaries/{sample}_boundaries.tsv",
        script = "bin/scripts/classify_deletions.py"
    output:
        inversion = "data/11_annotated_boundaries/{sample}_inversion.csv",
        deletion = "data/11_annotated_boundaries/{sample}_deletion.csv",
        # inversion_table = "data/11_annotated_boundaries/inversion_mechanism.csv",
        # deletion_table = "data/11_annotated_boundaries/deletion_mechanism.csv"
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
        inversion = "data/11_annotated_boundaries/{sample}_inversion.csv",
        script = "bin/scripts/inversion_replichore_classify.py"
    output:
        "data/11_annotated_boundaries/{sample}_inversion_classification.csv"
    params:
        input_dir = "data/11_annotated_boundaries/",
        ancestor = "Anc-_0gen_REL606",
        output_table = "inversion_replichores.csv"
    shell:
        """
        {input.script} --folder {params.input_dir} --oridif {input.ori_dif_coords} --ancestor {params.ancestor} --output {params.output_table}
        pwd
        """

# this rule just checks to see if the previous rule generated the main output tables. No shell executed

rule annotate_SV_mechanism_all:
    input:
        expand("data/11_annotated_boundaries/{sample}_inversion.csv", sample=df['assembly'].tolist()),
        expand("data/11_annotated_boundaries/{sample}_deletion.csv", sample=df['assembly'].tolist())
    output:
        inversion_table = "data/11_annotated_boundaries/inversion_mechanism.csv",
        deletion_table = "data/11_annotated_boundaries/deletion_mechanism.csv"

# god is watching you for this hideous code. repent for your sins.
        # """
        # cp {input.subject_path} {params.output_dir}/{params.subject_name} 
        # cp {input.query_path} {params.output_dir}/{wildcards.sample}.fasta
        # cd data/06_nucmer_alignment/{wildcards.sample}
        # touch ../../../{log}
        # nucmer --maxmatch -c 100 -b 500 -l 50 -p {wildcards.sample} {params.subject_name} {wildcards.sample}.fasta > ../../../{log} 2>&1
        # delta-filter -i {params.seq_id_cutoff} -l 100 {wildcards.sample}.delta > {wildcards.sample}.filtered.delta
        # show-coords -THrd {wildcards.sample}.filtered.delta > {wildcards.sample}.filtered.coords
        # touch {wildcards.sample}.done
        # rm {wildcards.sample}.fasta {params.subject_name}
        # cd ../../../
        # echo "Alignment complete. Working dir set to:"
        # pwd
        # """
