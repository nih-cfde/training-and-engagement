# Snakemake workflow tutorial

This tutorial is written for Unix or Linux compute environment (e.g., MacOS, Linux-based HPC).

## Introduction

### What is a workflow system and why use one?

Every computational workflow consists of multiple steps, starting with raw data and ending with summary data to plot and calculate statistics.

Workflow systems help you automate and manage the inputs, outputs, and commands for the analysis, thereby making it easier to maintain, reproduce, and share your workflow! Read more [here](https://www.nature.com/articles/d41586-019-02619-z).

### What is Snakemake?

[Snakemake](https://snakemake.readthedocs.io/en/stable/) is a Python-based workflow system ([see 2012 publication](https://academic.oup.com/bioinformatics/article/28/19/2520/290322)). The name 'Snakemake' comes from the fact that it's written in (and can be extended by) the Python programming language.

`snakemake` works by looking at a file, called a 'Snakefile', that contains rules for creating output files. Generally, each rule is defined as a step in the workflow. `snakemake` uses the rules and command line options to figure how the rules relate to each other so it can manage the workflow steps.

As an example, this tutorial will walk you through creating a `snakemake` workflow for variant calling. This tutorial was adapted from DIB lab course materials [here](https://github.com/ngs-docs/2020-GGG298) and [here](https://github.com/ngs-docs/2020-GGG201b-lab).

The contents of this tutorial are covered in three short videos:

- [Part 1](https://video.ucdavis.edu/media/snakemake+intro%2C+try+2/0_843yn8pn/166161802)

- [Part 2](https://video.ucdavis.edu/media/snakemake+intro+2+try+1/0_t1dpuzly)

- [Part 3](https://video.ucdavis.edu/media/snakemake+intro+3+try+1/0_gwnss4kq)

!!! info
    This may not be the variant calling workflow you would necessarily use in practice, but it serves as a good example for teaching Snakemake. Many people do indeed use `samtools`, but for particularly big or complex genomes, guidelines provided by  [GATK](https://gatk.broadinstitute.org/hc/en-us) would serve best. Additionally, various parameters associated with mapping, visualization etc may require tuning. 

The objectives of this tutorial are to:

!!! goal

    - learn how to write basic workflows with `snakemake` rules
    - learn variable substitution for `snakemake` rules
    - learn wildcard matching for `snakemake` rules
    - understand why workflow systems can help you do your computing more easily

## Set up

We will use conda to create an environment for this tutorial. If you don't have conda installed, please see the [Set up computing environment with conda on MacOS tutorial](../install_conda_tutorial.md).

!!! Tip

    Please refer to the [conda command cheatsheet](./conda_cheatsheet.md) for commonly used conda commands!

### 1. Download tutorial files:

First, create a new directory for this tutorial, e.g.,:
```
(base) $ mkdir learn_snakemake
```

We need two files for this tutorial. Click the links and save them in the directory you created above: 1) [environment.yml](./snakemake_tutorial_docs/environment.yml) and 2) [Snakefile](./snakemake_tutorial_docs/Snakefile.py).

Rename the `Snakefile.py` to `Snakefile`. There should be no file extension (we just added it so you'd be able to download the file!).

**TO DO: adding the .py extension enabled me to download the file by clicking the link, doesn't work when there is no extension. If you know how to do this better, please edit!**

### 2. Create new conda environment:

The environment.yml file tells conda 1) where to look for the software installations under 'channels' and 2) what software to install under 'dependencies'. You can also specify specific software versions, otherwise conda will download the most up-to-date version. Here are the specifications we'll use for this tutorial:
```
channels:
    - conda-forge
    - bioconda
    - defaults
dependencies:
    - bwa
    - snakemake-minimal=5.8.1
    - samtools=1.10
    - bcftools
```

Create environment:
```
(base) $ conda env create -n snaketest -f environment.yml
```

### 3. Activate conda environment:

```
(base) $ conda activate snaketest
```

### 4. Test that your environment is ready to go

You should have several software installed in your `snaketest` environment now. Check it out!
```
(snaketest) $ samtools --version
```
The output should say:
```
samtools 1.10
Using htslib 1.10.2
Copyright (C) 2019 Genome Research Ltd.
```

If you get an error, the software installation may have failed. You can check the software that is installed in your conda environment: `conda list`.

You can leave the conda environment with: `conda deactivate`.

Later in the tutorial, we'll use `wget` to download data. Installing `wget` on MacOS can be achieved by using `Homebrew`, a handy package installation manager. This step will take a few minutes and the installation is outside the conda environment `snaketest`:

```
# install Homebrew
(base) $ /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

# use brew command to install wget
(base) $ brew install wget

# test installation
(base) $ wget --version
```
Alternatively, we can also use conda to install `wget`, which will work in both Unix and Linux based systems:

```
conda install -c anaconda wget
```
