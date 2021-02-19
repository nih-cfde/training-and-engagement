---
layout: page
title: RNAseq Tutorial Overview
---

Differential Gene Expression Analysis on Cavatica Cloud Platform
====================================================================

**RNA sequencing (RNAseq)** is a high throughput technique that provides qualitative and quantitative information about RNA biology including transcriptome-wide expression quantification, discovery of novel genes and gene isoforms, and differential expression.

The goal of this tutorial is to enable you to: </br>

1. create virtual cancer cohorts using the NIH Common Fund-supported Gabriella Miller Kids First Data portal (KF Portal). </br>
2. analyze the differential gene expression (DGE) on Cavatica, an integrated cloud based platform.

You will learn two different approaches for DGE analysis using open access human cancer data on Cavatica: **(a)** using a public workflow app
**(b)** running code from an analysis script on an instance with RStudio computational environment.

**Table of contents**

| Est. Time| Lesson Name | Description|
| ---|--------|--------|
| 10 mins |[An Introduction to RNA-Seq](./rna_seq_2.md)| Background about RNAseq
| 20 mins |[Selecting Kids First Cancer Cohort](./rna_seq_3.md)| Select Kids First open access cancer RNAseq files and push to Cavatica  |
| 20 mins |[Cavatica - View, Filter, Tag and Download](./rna_seq_4.md) | Filter imported data, tag and download relevant metadata from Cavatica |
| 20 mins |[Setup DESeq2 Public App](./rna_seq_5.md)| Setting up the workflow app based on DESeq2 on Cavatica |
| 15 mins |[Phenotype File and Upload to Cavatica](./rna_seq_6.md) | Reformat metadata file and upload it to Cavatica |
| 50 mins |[Analysis with DESeq2 Public App](./rna_seq_7.md) | Run the DESeq2 app with appropriate inputs and computational settings |
| 60 mins |[Analysis using Data Cruncher](./rna_seq_8.md) | Analysis on an instance in the RStudio environment |

!!! note "Learning Objectives"
    * learn to build virtual cohorts on KF portal
    * learn to navigate project folder and perform file operations on Cavatica
    * learn to upload and download data from Cavatica
    * learn to search, copy, and edit public workflow apps on Cavatica
    * learn to perform differential gene expression (DGE) analysis using DESeq2 app
    * learn to setup analysis environment and execute code for DGE analysis

=== "Prerequisites"
    * Setup: Integrated login accounts on Kid's First Data Portal & Cavatica - Follow our lessons on account setup and connecting the two accounts.
           - [Setup Kids First account](../Kids-First/Portal-Setup-And-Permissions/KF_3_KF_Registration.md)
           - [Register for Cavatica](../Kids-First/Portal-Setup-And-Permissions/KF_4_Cavatica_Registration.md)
           - [Connect KF and Cavatica](../Kids-First/Portal-Setup-And-Permissions/KF_5_ConnectingAccounts.md)

    !!! note "Login Credentials"

        You do not need **eRA Commons ID** to do the lesson!

    * Background: Knowledge of biology and rudimentary genetics.
    * Technology: Basic knowledge of R and command line. Familiarity with RStudio is useful.
    * Financial: Pilot funds ($100) are provided to every user on Cavatica with linked KF accounts.
    * Time: Initial account setup may take hours to a day for verification. Setup of eRA Commons ID may take days and is institute dependent.
=== "Est.Cost"
    - DESeq2 app < $1.00
    - Analysis with R < $1.00      
=== "Tutorial Resources"
    - [Kids First Data Portal](https://kidsfirstdrc.org)
    - [Cavatica Documentation](https://docs.cavatica.org/docs/getting-started)
    - [Playlist of video tutorials explaining concepts used in RNAseq analysis](https://www.youtube.com/playlist?list=PLblh5JKOoLUJo2Q6xK4tZElbIvAACEykp)
    - [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#how-do-i-use-vst-or-rlog-data-for-differential-testing)
    - [tximport](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)
