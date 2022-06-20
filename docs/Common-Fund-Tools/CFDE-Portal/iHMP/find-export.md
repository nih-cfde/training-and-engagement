---
layout: page
title: Using the CFDE Search Portal to Find Files
---

#  Using the CFDE Search Portal to Find Files

The Common Fund Data Ecosystem Coordinating Center (CFDE-CC) supports efforts to make Common Fund data sets more findable, accessible, interoperable, and reusable (FAIR) for the scientific community through collaboration, end-user training, and data set sustainability. The CFDE-CC manages and organizes CFDE activities, engages with participating Common Fund programs, connect with user communities, supports training, develops tools and standards, and provides technical expertise to Common Fund programs.

The [CFDE Search Portal](https://app.nih-cfde.org/) uses the Crosscut Metadata Model (C2M2) a flexible metadata standard for describing experimental resources in biomedicine and related fields. This portal supports faceted search of metadata concepts such as anatomical location, species, and assay type, across a wide variety of datasets using a controlled vocabulary (we do not currently support protected metadata). This allows researchers to find a wide variety of data that would otherwise need to be searched individually, using varying nomenclatures. The portal only accepts C2M2 data packages from Common Fund Programs.

This tutorial focuses on the Human Microbiome Project (HMP). The goal of this tutorial is to identify small FASTQ files with persistent identifiers from a longitudinal multi-comics study. Please refer to the [Portal User Guide](https://docs.nih-cfde.org/en/latest/about/portalguide/) for a detailed description of all portal features. 

!!! note "Learning Objectives"
    
    -   learn how to access the CEFD search portal
    -   learn how to create personal collections
    -   learn how to search for files meeting a specific criterion

## Create a new personal collection

Go to [the CFDE data portal at app.nih-cfde.org/](https://app.nih-cfde.org/). 

Log in (upper right).

![](https://i.imgur.com/YVXgMVK.png)

Under your username (upper right), create a new personal collection. 

![](https://i.imgur.com/D2eEXg2.png)

For name, you can use "Tuesday demo" or anything else. You can leave description blank.

![](https://i.imgur.com/eNoJFep.png)


### Find some files

Go back to [the CFDE data portal main page](https://app.nih-cfde.org/). 

Select <span class="highlight_txt">File</span> (upper left).

![](https://i.imgur.com/nIlZ2Jw.png)

Use the facets on the left to select:
* Common Fund Program: HMP
* Project: "Longitudinal multi 'omics"
* "has persistent ID" - True
* Uncompressed size in Bytes - 50000000 to 60000000 (50 MB to 60 MB).

![](https://i.imgur.com/7SAZK0X.png)


![](https://i.imgur.com/9wOPGAY.png)


With these selections, the first result should have "Filename" of `SRR5935743_1.fastq`.

![](https://i.imgur.com/ULbqD7W.png)


## Add files to your personal collection

Click on the first result to get a detail view. Then add it to your personal collection:
* scroll down to "Part of personal collection" and click  <span class="highlight_txt">Link records</span>.
* Select your personal collection, click <span class="highlight_txt">Link</span> (upper right).

![](https://i.imgur.com/76nE4vc.png)

![](https://i.imgur.com/lWWBz0m.png)


Click <span class="highlight_txt">back</span> in your browser, to get back to your filtered search.

Repeat linking to a collection with the third result (Filename: `SRR5950647_1.fastq`). Add it to the same personal collection.

![](https://i.imgur.com/aOX8Gz7.png)


!!! info
    For today, here are the direct links to the two files we'll be using:    
    -   [file 1, `SRR5935743_1.fastq`](https://app.nih-cfde.org/chaise/record/#1/CFDE:file/nid=528892)
    -   [file 2, `SRR5950647_1.fastq`](https://app.nih-cfde.org/chaise/record/#1/CFDE:file/nid=531342)


(Note, you could select any files you like, but these are small enough to work and I know what the results will be. So it's good for today's demo; I suggest trying new/different files as a Thursday exercise!)

## Export your personal collection

Go to your collection, and select <span class="highlight_txt">export</span> and choose <span class="highlight_txt">NCPI manifest format</span>.

!!! info
    **What is NCPI?**
    "NCPI" stands for "NIH Cloud Platform Interoperability", an effort by the NIH to convene around interoperation for cloud workbenches.
 

### Examine the NCPI manifest file

You should now have a CSV file in your Downloads that, when examined with a spreadsheet program, looks like this:

![](https://hackmd.io/_uploads/HkshxgdLq.png)

The key piece of information in here is the `drs_uri` column, which provides a [Data Repository Service](https://ga4gh.github.io/data-repository-service-schemas/preview/release/drs-1.0.0/docs/) location from which to download the files.
