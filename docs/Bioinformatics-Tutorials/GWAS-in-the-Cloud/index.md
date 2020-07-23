---
layout: page
title: GWAS tutorial overview
---

GWAS in the cloud using Amazon Web Services
===========================================

This is a detailed tutorial on performing genome wide association analysis (GWAS) using Amazon Web Services (AWS) based on the [ANGUS 2017 GWAS tutorial](https://angus.readthedocs.io/en/2017/GWAS.html)

For this tutorial, we will **not** work with human data. We will use coat color in dogs as the trait of interest (instead of disease), and test the association of a genome-wide set of SNPs with two coat color variants: yellow and dark. To extrapolate this tutorial to human disease data, you might consider yellow coat color phenotype as the "case" (or disease) and dark coat color as the "control" (or normal) condition.

The goal of this tutorial is to look at summary statistics and perform a simple association analysis starting with a variant calling file (.vcf). We will also produce Manhattan plots to visualize loci associated with yellow coat coloration in dogs.


**Table of contents**

- [What is GWAS?](background.md)
- [Tutorial Overview](tutorial_overview.md)
- [Setting up an AWS instance](aws_instance_setup.md)
- [Accessing the AWS instance](Accessing_aws.md)
- [Install PLINK](plink_install.md)
- [Install VCFtools](vcftools_install.md)
- [Install R and RStudio](RStudio.md)
- [Download the data](data_download.md)
- [Analyze](analyze.md)
- [Manhattan Plots](manhattan.md)
- [Terminate AWS Instance](terminate_aws.md)
