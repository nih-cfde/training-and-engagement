---
layout: page
title: GWAS tutorial overview
---

How to do GWAS in the cloud using Amazon Web Services
=====================================================

**Genome-wide association studies (GWAS)** offer a way to rapidly scan entire genomes and find genetic variation associated with a particular disease condition.

Our aim is to teach researchers how to perform genome wide association analysis using Amazon Web Services (AWS). This tutorial will enable researchers with minimal bioinformatics background to set up and access an AWS instance, move data in and out of the AWS instance, run some basic summary statistics and perform a simple association analysis starting with variant calling files (.vcf). We will also produce Manhattan plots to visualize variants associated with traits.

For this tutorial, we will *not* work with human data. We will use coat color in dogs as the trait of interest (instead of disease), and test the association of a genome-wide set of single nucleotide polymorphisms (SNPs) with two coat color variants: yellow and dark. To extrapolate this tutorial to human disease data, you might consider yellow coat color phenotype as the "case" (or disease) and dark coat color as the "control" (or normal) condition.

This tutorial is based on the [ANGUS 2017 GWAS tutorial](https://angus.readthedocs.io/en/2017/GWAS.html)

**Table of contents**

| Est. Time| Lesson Name | Description|
| ---|--------|--------|
| 10 mins |[What is GWAS?](background.md)| Background                   
| 20 mins |[Set up an AWS instance](aws_instance_setup.md)|How to set up an amazon web services instance|
| 40 mins |[Download and move data to AWS](download_accessAWS.md) | Download dog coat color data and use terminal to access data on remote computer |
| 10 mins |[Install PLINK](plink_install.md)| Install the software PLINK |
| 10 mins |[Install VCFtools](vcftools_install.md) | Install the software vcftools |
| 10 mins |[Install R and RStudio](RStudio.md) | Install the software R and RStudio |
| 10 mins |[Analyze](analyze.md) | Generate summary statistics and association analysis |
| 20 mins |[Manhattan Plots](manhattan.md) | Make some plots to visualize data |
| 10 mins |[Terminate AWS Instance](terminate_aws.md) | Shut down the cloud computer |

!!! note "Learning Objectives"
    * learn to set up and access Amazon Web Services
    * learn to move data in and out of the AWS instance
    * learn to install and run all software necessary for GWAS analysis
    * learn to produce Manhattan plots

=== "Prerequisites"
    * Background: Some expertise in biology and fundamental genetics.
    * Technology: Basic shell scripting knowledge and access to MacOS. Users must be comfortable with finding and opening a terminal window, navigating to specific directories and running pre-scripted commands in the terminal.
    * Financial: First time AWS users require a valid credit card to set up an AWS account.
    * Time: AWS account setup needs approval by AWS, and approval times can range from minutes to days.
=== "Tutorial Resources"
    - [AWS website](http://aws.amazon.com/)

    - [Manhattan plots explained](https://www.google.com/search?q=how+to+read+a+manhattan+plot&oq=how+to+read+a+manhattan+plot&aqs=chrome..69i57.7911j0j4&sourceid=chrome&ie=UTF-8#kpvalbx=_tXIPX9mmFsmT0PEP64-OkAk26)
