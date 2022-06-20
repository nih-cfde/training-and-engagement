---
layout: page
title: Movement Disorders Datasets in CF Program Portals
---

<div class="banner"><span class="banner-text">This tutorial is outdated. Please refer to the Common Fund Tools page for up-to-date tutorials.</span></div>

# Movement Disorders Datasets in CF Program Portals

[Using the `csv` manifest details we just exported from the CFDE portal](./movement-disorders-portal-export.md), the associated data files can be accessed from data portals of the individual Common Fund programs, [Metabolomics](https://www.metabolomicsworkbench.org) and [LINCS](http://lincsportal.ccs.miami.edu/datasets/) respectively. In this use case, the selected cohort from the CFDE Portal represents genomic, transcriptomic and metabolomic data from induced pluripotent stem cells and neural stem cells of individuals with ALS, SMA (LINCS), Parkinson's and MS (Metabolomics) along with healthy controls.


## Metabolomics WorkBench

The details of each Metabolomics study can be viewed, analyzed and downloaded using the [Metabolomics WorkBench](https://www.metabolomicsworkbench.org). The `persistant_id` for each study is associated with a summary page which lists all the available analyzed and raw data, metadata associated with study design, experimental conditions, sample preparation details and analysis techniques along with   contributor information and creation date.

![Metabolomics Workbench summary](../images/Metabolomics-workbench-summary.png "Metabolomics Workbench summary")

The associated metadata for the different fields can be listed as tabs. Study data can be downloaded as `zip` files using the "Download data" option. Some studies have processed datasets while some have raw and processed datasets. One can use the "View archive contents" to explore the contents of the `zip` file prior to download.

![Metabolomics data download](../images/Metabolomics-download-data.png "Metabolomics data download")

Selecting the `Perform statistical analysis` lists multiple options for statistical tests, clustering, pathway mapping and visualization that can be run on the study data.

![Metabolomics statistical tests](../images/Metabolomics-statistical-test.png "Metabolomics statistical tests")

## LINCS Data Portal

The `persistant_id` for each LINCS dataset is linked to the study page in the [LINCS data portal](http://lincsportal.ccs.miami.edu/datasets/) which list the Description, Metadata and Download tabs. The associated metadata and analyzed data are available for direct download.

![LINCS data portal](../images/LINCS-data-portal.png "LINCS data portal")

The "Data Source" lists the link to the dbGaP study which all relevant study description, sequencing details, associated published literature along with information to apply for access to the controlled experimental data.

![LINCS data source](../images/LINCS_dbGaP-data-source.png "LINCS data source")
