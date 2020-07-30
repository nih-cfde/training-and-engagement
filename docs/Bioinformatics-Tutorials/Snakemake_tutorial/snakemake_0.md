# Snakemake Introduction

## What is a workflow system and why use one?

Every computational workflow consists of multiple steps, starting with raw data and ending with summary data to plot and calculate statistics.

Workflow systems help you automate and manage the inputs, outputs, and commands for the analysis, thereby making it easier to maintain, reproduce, and share your workflow! Read more [here](https://www.nature.com/articles/d41586-019-02619-z).

## What is Snakemake?

[Snakemake](https://snakemake.readthedocs.io/en/stable/) is a Python-based workflow system ([see 2012 publication](https://academic.oup.com/bioinformatics/article/28/19/2520/290322)). The name 'Snakemake' comes from the fact that it's written in (and can be extended by) the Python programming language.

Snakemake works by looking at a file, called a 'Snakefile', that contains rules for creating output files. Generally, each rule is defined as a step in the workflow. Snakemake uses the rules and command line options to figure how the rules relate to each other so it can manage the workflow steps.

Let's get started!
