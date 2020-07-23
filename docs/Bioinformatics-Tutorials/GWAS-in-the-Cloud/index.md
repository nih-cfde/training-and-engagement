---
layout: page
title: GWAS tutorial overview
---

GWAS in the cloud using Amazon Web Services
===========================================

Genome-wide association studies (GWAS) offer a way to rapidly scan entire genomes and find genetic variation associated with a particular disease condition.

Our aim is to teach researchers how to performing genome wide association analysis using Amazon Web Services. This tutorial will enable researchers with minimal bioinformatics background to set up and access an AWS instance, move data in and out of the AWS instance, run some basic summary statistics and perform a simple association analysis starting with variant calling files (.vcf). We will also produce Manhattan plots to visualize variants associated with traits.

For this tutorial, we will **not** work with human data. We will use coat color in dogs as the trait of interest (instead of disease), and test the association of a genome-wide set of single nucleotide polymorphisms (SNPs) with two coat color variants: yellow and dark. To extrapolate this tutorial to human disease data, you might consider yellow coat color phenotype as the "case" (or disease) and dark coat color as the "control" (or normal) condition.

This tutorial is based on the [ANGUS 2017 GWAS tutorial](https://angus.readthedocs.io/en/2017/GWAS.html)



**Table of contents**


| Est. time| Lesson name | Description|
| ---|--------|--------|
| 00:05|[What is GWAS?](background.md)| Background                   
| 00:15|[Set up an AWS instance](aws_instance_setup.md)|How to set up an amazon web services instance|
| 00:15| [Access the AWS instance](Accessing_aws.md) | Use terminal to connect to the remote computer |
| 00:05| [Install PLINK](plink_install.md)| Install the software PLINK |
| 00:05| [Install VCFtools](vcftools_install.md) | Install the software vcftools |
| 00:05| [Install R and RStudio](RStudio.md) | Install the software R and RStudio |
| 00:15| [Download the data](data_download.md) | Download dog coat color data into AWS instance |
| 00:30| [Analyze](analyze.md) | Generate summary statistics and association analysis |
| 00:20| [Manhattan Plots](manhattan.md) | Make some plots to visualize data |
| 00:02| [Terminate AWS Instance](terminate_aws.md) | Shut down the cloud computer |
