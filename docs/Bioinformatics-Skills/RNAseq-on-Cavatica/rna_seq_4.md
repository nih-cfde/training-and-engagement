---
layout: page
title: Cavatica - View, Filter, Tag and Download
---

To view the project folder on Cavatica, you can click the link from the pop box in KF portal after successful copy of files which will open the Cavatica login page.
Alternatively, you can [login to Cavatica](https://cavatica.sbgenomics.com){:target="_blank"} in a new tab.

## Step 1: View files in Cavatica <a name="view-files"></a>

* Select the newly created project folder under the <span class="highlight_txt">Projects</span> tab.
* The Dashboard of the project folder has three panels: Description, Members and Analyses.
* Click on the <span class="highlight_txt">Files</span> tab to list all the project files.

![Files tab in project homepage](../rna-seq-images/10_Cavatica.png "Files tab in project homepage")

* Click on the <span class="highlight_txt">Type: All</span> filter for a drop down box which lists the type and number of files: 99 compressed tsv files.

![Total files in project folder](../rna-seq-images/11_Cavatica.png "Total files in project folder")


## Step 2: Apply filters to subset cohort

Before we proceed to the Differential Gene Expression Analysis (DGE analysis), it is a good idea to examine the metadata associated with our selected cohort. Since we aim to keep the experimental design simple, we will further filter down to remove possible sources of variation.

The columns visible in the table are the platform default options. Click on <span class="highlight_txt">:fontawesome-solid-columns:</span> on the right hand corner and select any columns to view from the metadata list.

![Edit table columns](../rna-seq-images/12_Cavatica.png "Edit table columns")

![Custom columns](../rna-seq-images/13_Cavatica.png "Custom columns"){: align=right width=52%}

Here we have selected:

 - Age at diagnosis
 - Vital status
 - tumor_location
 - histology
 - histology_type

</br>
!!! info "Age at diagnosis"

    The default unit for any age metadata field is recorded in days and is reflected in the large numeric values for Age at diagnosis column.

Each of these columns have multiple values. To filter the data using values within multiple metadata columns, use the <span class="highlight_txt">:fontawesome-solid-plus:</span> sign to add a filter. If you cannot see the <span class="highlight_txt">:fontawesome-solid-plus:</span> button, refresh your browser, as your session may have timed out.

![Apply additional filters](../rna-seq-images/14_Cavatica.png "Apply additional filters")

* First, we filter to only include surviving patients. Click on <span class="highlight_txt">:fontawesome-solid-plus:</span> and
choose <span class="highlight_txt">Vital status</span>, then select **Alive** from the sub-menu.

![Vital status filter](../rna-seq-images/15_Cavatica.png "Vital status filter")

* Since the patients could have presented with multiple cancers over diagnostic timeline, the <span class="highlight_txt">histology</span> metadata has other values in addition to the cancer types of interest. Click <span class="highlight_txt">:fontawesome-solid-plus:</span> again this time choosing <span class="highlight_txt">histology</span> and selecting both **Medulloblastoma** & **Ependymoma**.

![histology filter](../rna-seq-images/16_Cavatica.png "histology filter")

* To ensure comparison of cancer from the first presentation in the patient, we eliminate recurrent or progressive subtypes using the <span class="highlight_txt">histology_type</span> filter following the same steps as previously. This time select only **Initial CNS Tumor**.

![histology_type filter](../rna-seq-images/17_Cavatica.png "histology_type filter")

The tumor_location metadata column has some values that include multiple anatomically distinct locations separated by a **`;`**. This could indicate the observation of spread of tumor to multiple locations during first occurrence.

* We filter using the <span class="highlight_txt">tumor_location</span> metadata, choosing only values **without** the **`;`**. Select the eleven distinct values for tumor_location (not including those with **`;`**, **`Not Reported `** , and **`Other locations NOS `**). You can see the complete list in the screen capture below.

![tumor_location filter](../rna-seq-images/18_Cavatica.png "tumor_location filter")

This results in total of 50 files from our initial 98 copied files.

## Step 3: Create tags & download filtered dataset

To enable quick access to the filtered data without having to re-run all the metadata filters, we can create tags for the filtered data.


* Select all the files by clicking on <span class="highlight_txt">:material-square-rounded-outline:</span> in the column header and click on <span class="highlight_txt">:fontawesome-solid-tags:Tags</span> tab.  

![All filtered files](../rna-seq-images/19_Cavatica.png "All filtered files")

* Type the name of the tag and click <span class="highlight_txt">Add new tag</span>.

![Add new tag](../rna-seq-images/20_Cavatica.png "Add new tag")

!!! tip "Tag Names"

    While you can use any tag name you see fit, use **DGE-FILTER-DATA** as used in this lesson to match your screen with the lesson screenshots.

* Click <span class="highlight_txt">Apply</span>. In case, you wish to remove the tag, use the <span class="highlight_txt">:fontawesome-solid-times:</span> in the tag name to delete.

![Apply new tag](../rna-seq-images/21_Cavatica.png "Apply new tag")

The filtered files are now tagged. We need to download and modify the metadata file which will be used as the accompanying phenotype file for our DGE analysis in the next lesson. To download:

* Click on the <span class="highlight_txt">:fontawesome-solid-ellipsis-h:</span> button on the right corner.
* Select <span class="highlight_txt">Export metadata manifest from filtered files</span>.

![Download filtered metadata](../rna-seq-images/22_Cavatica.png "Download filtered metadata")

In our next lesson, we will learn to setup the DESeq2 app in our project folder.
