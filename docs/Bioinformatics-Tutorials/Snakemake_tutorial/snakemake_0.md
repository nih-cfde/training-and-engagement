# Snakemake Introduction

This tutorial is written for Unix or Linux compute environment (e.g., MacOS, Linux-based HPC).

## Introduction

### What is a workflow system and why use one?

Every computational workflow consists of multiple steps, starting with raw data and ending with summary data to plot and calculate statistics.

Workflow systems help you automate and manage the inputs, outputs, and commands for the analysis, thereby making it easier to maintain, reproduce, and share your workflow! Read more [here](https://www.nature.com/articles/d41586-019-02619-z).

### What is Snakemake?

[Snakemake](https://snakemake.readthedocs.io/en/stable/) is a Python-based workflow system ([see 2012 publication](https://academic.oup.com/bioinformatics/article/28/19/2520/290322)). The name 'Snakemake' comes from the fact that it's written in (and can be extended by) the Python programming language.

Snakemake works by looking at a file, called a 'Snakefile', that contains rules for creating output files. Generally, each rule is defined as a step in the workflow. Snakemake uses the rules and command line options to figure how the rules relate to each other so it can manage the workflow steps.

As an example, this tutorial will walk you through creating a Snakemake workflow for variant calling. This tutorial was adapted from DIB lab course materials [here](https://github.com/ngs-docs/2020-GGG298) and [here](https://github.com/ngs-docs/2020-GGG201b-lab).

The contents of this tutorial are covered in three short videos:

- Part 1: [Introducing the Snakefile and Snakemake](./snakemake_2.md)

- Part 2: [Continue decorating the Snakefile](./snakemake_4.md)

- Part 3: [Running through the entire Snakemake workflow](./snakemake_4.md)

!!! info
    This may not be the variant calling workflow you would necessarily use in practice, but it serves as a good example for teaching Snakemake. Many people do indeed use `samtools`, but for particularly big or complex genomes, guidelines provided by  [GATK](https://gatk.broadinstitute.org/hc/en-us) would serve best. Additionally, various parameters associated with mapping, visualization etc may require tuning.

The objectives of this tutorial are to:

!!! goal

    - learn how to write basic workflows with Snakemake rules
    - learn variable substitution for Snakemake rules
    - learn wildcard matching for Snakemake rules
    - understand why workflow systems can help you do your computing more easily
