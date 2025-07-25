---
title: 'seabreeze: A Pipeline for Analyzing Structural Variation Between Bacterial Genome Assemblies'
tags:
  - structural variation
  - experimental evolution
  - bacterial genomics
  - Snakemake
authors:
  - name: Ira Zibbu
    orcid: 0009-0003-3930-5397
    affiliation: "1, 2"
  - name: Claus O. Wilke
    orcid: 0000-0002-7470-9261
    affiliation: 3
  - name: Jeffrey E. Barrick
    orcid: 0000-0003-0888-7358
    corresponding: true
    affiliation: 1
affiliations:
  - name: Department of Molecular Biosciences, The University of Texas at Austin, Austin, Texas, United States of America
    index: 1
  - name: School of Biology, Indian Institute of Science Education and Research, Thiruvananthapuram, Kerala, India
    index: 2
  - name: Department of Integrative Biology, The University of Texas at Austin, Austin, Texas, United States of America
    index: 3
date: 31 January 2025
bibliography: paper.bib
---

# Summary

Structural mutations—such as large insertions, deletions, duplications, inversions, and translocations—play a unique and important role in bacterial evolution. Recent advances in long-read sequencing have made accurate, high-throughput predictions of structural mutations possible. *seabreeze* is a tool for comprehensively analyzing genetic variation among bacterial genomes caused by structural mutations. It manages a workflow that combines existing packages and custom scripts to automate and unite several analyses into a single, easy-to-use pipeline. For specified pairs of bacterial genomes, *seabreeze* predicts and visualizes structural mutations, and annotates the affected genes. It also uses information about transposons at the boundaries of structural mutations to predict their involvement in generating the mutations. Finally, *seabreeze* provides information about size differences between the two genomes and changes in replichore balance within circular chromosomes. Although *seabreeze* was developed to characterize structural mutations that evolved in laboratory experiments, it can be used to analyze any sufficiently closely related genomes from strains of the same bacterial species.

# Statement of Need

Evolution experiments are one of many approaches that are used to understand the processes that shape microbial genomes. In these experiments, populations of microbes are propagated under controlled conditions for a sufficient period of time to observe evolution in action. Then, evolved genomes are resequenced and compared to the ancestral genome. We created *seabreeze*, a Snakemake pipeline that automates comparing genome assemblies to predict and analyze structural variation that emerges in these experiments. *seabreeze* fulfills two core needs. First, although a wide range of software programs exist to identify structural variation from high-throughput sequencing data [@ahsan_survey_2023], they make assumptions that are only appropriate for eukaryotic genomes and/or are limited by how they compare reads to a reference genome. In contrast, *seabreeze* is explicitly tailored for bacterial genome analysis and takes advantage of the benefits of comparing genome assemblies. Second, *seabreeze* unites several standalone open-source tools and new custom scripts into a single easy-to-use pipeline for comprehensive bacterial genome analysis. Other notable tools for bacterial genome analysis such as Artemis [@carver_act_2005] and Mauve [@darling_progressivemauve_2010] can detect and visualize rearrangements in bacterial genomes but lack other functionality, such as predicting the mechanisms of structural mutations and annotating what genes they affect.

# Implementation

The latest release of *seabreeze* can be downloaded from the [official GitHub repository](https://github.com/barricklab/seabreeze). The [online documentation](https://barricklab.github.io/seabreeze/) details installation, usage, and output. It contains a tutorial to demonstrate how it can be used. *seabreeze* uses Snakemake [@koster_snakemakescalable_2012] to manage the pipeline and automatically install and manage dependencies for the workflow from `conda-forge` [@conda-forge_community_conda-forge_2015] and `bioconda` repositories [@the_bioconda_team_bioconda_2018].  Users supply fully assembled sequences in FASTA format for each reference-query pair, and specify which pairs of sequences to compare in a CSV file to allow for batch-processing. Suitable input files can be generated from long-read DNA sequencing datasets with assemblers such as Flye [@kolmogorov_assembly_2019], Canu [@koren_canu_2017] or Raven [@vaser_time-_2021], or from consensus assembly tools like Trycycler [@wick_trycycler_2021]. 

*seabreeze* has several subcommands to perform specific tasks. It begins by accounting for the circular nature of bacterial chromosomes and plasmids by rotating the reference-query pairs to a common start coordinate. Next, it uses the following packages to perform various analysis steps for each reference-query pair: (1) compute size difference between genomes with `biopython` [@cock_biopython_2009]; (2) predict the locations of transposons with `ISEScan` [@xie_isescan_2017]; (3) align genomes with `mummer4` [@marcais_mummer4_2018]; (4) predict structural mutations with `SyRI` [@goel_syri_2019] and filter false-positive calls with a custom script; (5) generate intuitive synteny plots to visualize structural variation with a customized version of `plotsr` (Figure 1) [@goel_plotsr_2022]; (6) annotate the genes contained in the mutated regions with `prokka` [@seemann_prokka_2014] and `breseq` [@sun_identification_2014]. 

![Synteny plot generated by seabreeze. This plot compares an ancestor genome (top, in blue) to its simulated evolved genome (bottom, in orange). Grey regions are syntenic (i.e., gene presence and order is preserved between both genomes). A single large inversion (orange ribbon) and several deletions (red ribbons) are visible.\label{fig1}](REL606_evolved_1.png)

Most bacterial chromosomes are circular and contain a single origin and terminus of replication, which are used to define an origin-terminus axis that divides the genome into two halves called replichores [@rocha_replication-related_2004]. *seabreeze* introduces new scripts to analyze the placement of inversions relative to the origin-terminus axis and how they affect the symmetry of the two replichores. Most structural mutations occur through recombination between homologous sequences, and in particular, bacterial genomes often contain multiple copies of simple transposons (also known as insertion sequences) that serve as sites for recombination [@achaz_associations_2003]. *seabreeze* uses the locations of transposons to predict whether they were involved in generating inversions and deletions. Highly diverged genomes may have successive structural mutations or low sequence homology in syntenic regions that complicate these inferences, and *seabreeze* may be unable to identify the individual mutational steps that led to complex structural variation in these cases.

Overall, the growing use of long-read sequencing will increasingly make it possible to  compare evolved bacterial genomes to their ancestors at the assembly level versus by mapping reads to a reference genome. *seabreeze* will be of utility to researchers investigating how and why the organization of bacterial genomes evolves.


# Acknowledgments 

Development of *seabreeze* was supported by the National Science Foundation (DEB-1951307). We acknowledge the Texas Advanced Computing Center (TACC) at the University of Texas at Austin for providing computing resources. Ira Zibbu acknowledges additional support from the Khorana Scholars Program and the DST INSPIRE Fellowship. Claus Wilke acknowledges support from the Blumberg Centennial Professorship in Molecular Evolution at the University of Texas at Austin. We thank the developers of `plotsr` for making their code available under the MIT license.

# References
