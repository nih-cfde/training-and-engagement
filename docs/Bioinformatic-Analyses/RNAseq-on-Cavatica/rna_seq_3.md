---
layout: page
title: Selecting a Kids First Cancer Cohort
---

<div class="banner"><span class="banner-text">Lesson in Development</span></div>

Selecting a Kids First Cancer Cohort
====================================

The [Gabriella Miller Kids First Pediatric Data Portal (KF Portal)](https://kidsfirstdrc.org) hosts datasets at the intersection of congenital birth defects and pediatric cancers, with genomic files for more than 16,000 participants.

!!! tip "Kids First Data Portal"

    Check out our lessons on Kids First to learn more about the [different Data Portal features](../../Common-Fund-Tools/Kids-First/Exploring-Data-in-the-KF-Portal/KF_5_Explore.md) and [building simple to complex queries](../../Common-Fund-Tools/Kids-First/Advanced-KF-Portal-Queries/KF_9_AdvancedQuery.md).

Files on the KF Portal are managed through different access levels. _Open access_ files (including processed files of somatic samples) can be viewed and downloaded by any user. _Controlled access_ files (including raw sequencing files and imaging data) require approvals through dbGaP. For this tutorial, we will use **open access pre-processed files** generated using [Kallisto (v0.43.1)](http://pachterlab.github.io/kallisto//releases/2017/03/20/v0.43.1), which uses pseudoalignments to quantify transcript abundance from raw data.

!!! info "KFDRC RNAseq Workflow"

    The [Kids First RNAseq Workflow](https://github.com/kids-first/kf-rnaseq-workflow) uses multiple tools/packages for expression detection and fusion calls. The workflow requires raw FASTQ files (controlled access) as input and generates multiple outputs including the Kallisto transcript quantification files. All the output files of this pipeline are available on the KF Portal as open access data. Due to access restrictions and computational intensity, this tutorial will not cover the Kallisto workflow, but users with their own RNAseq data may consider starting from this point.

## Step 1: Filter for open access data

* Login to the [KF Portal](https://kidsfirstdrc.org/).
* Select the <span class="highlight_txt">File Repository</span> tab.
* Select the <span class="highlight_txt">Browse All</span> option for the Filter.

![File Repository](../rna-seq-images/1_KFDRC.png "File Repository")

!!! note "Data Summary"

    At the time this tutorial was written (January 2021), the portal contained a total of 88,728 files. Because new datasets are constantly uploaded to the KF Portal, exact numbers within your query may change slighly when run in the future.

* Select the <span class="highlight_txt">Access</span> filter listed under <span class="highlight_txt">FILE</span> field.
* Select <span class="highlight_txt">Open</span> value.
* Click <span class="highlight_txt">View Results</span> to update selection. This results in 18,162 files.

![Open access filter](../rna-seq-images/2_KFDRC.png "Open access filter")

## Step 2: Apply File Filters to obtain RNAseq files

Select the <span class="highlight_txt">File Filters</span> tab and apply the following filters:

* **Experimental Strategy**: `RNA-Seq`
* **Data Type**: `Gene expression`
* **File Format**: `tsv`

This results in 1,477 files.

![File filters](../rna-seq-images/3_KFDRC.png "File filters")

## Step 3: Select cancer type

Switch to the <span class="highlight_txt">Clinical Filter</span> tab and apply the following.

* **Diagnosis (Source Text)**: `Medulloblastoma` and `Ependymoma`

This filters the number of files to 235.

![Cancer type](../rna-seq-images/4_KFDRC.png "Cancer type")

## Step 4: Subset cohort

To reduce possible sources of variation due to participant demographics, we will further narrow the query to only include data from white male patients.

Under the <span class="highlight_txt">Clinical Filters</span> tab select:

* **Gender**: `Male`
* **Race**: `White`

This results in 99 files.

![Subset by Clinical Filters](../rna-seq-images/5_KFDRC.png "Subset by Clinical Filters")

## Step 5: Copy files to Cavatica

!!! important "Important"

    It is crucial to ensure the Cavatica integrations are enabled to allow for file transfers. Find more details in our [Push to Cavatica lesson](../../Common-Fund-Tools/Kids-First/KF_7_PushToCavatica.md). You **do not** have to have the Data Repository Integrations set up to continue with this lesson.

* Click on the <span class="highlight_txt">ANALYZE IN CAVATICA</span> button.
* Select the <span class="highlight_txt">CREATE A PROJECT</span> option and provide an appropriate name for your folder. This tutorial uses **`cancer-dge`** as the project name.
* Use the <span class="highlight_txt">SAVE</span> option to create the project.

![Create project on Cavatica](../rna-seq-images/6_KFDRC.png "Create project on Cavatica")

Following project creation, the option will update to enable copying of the selected files to Cavatica.

![Copy files to Cavatica](../rna-seq-images/7_KFDRC.png "Copy files to Cavatica")

Successful copying of the files to the project folder will result in a pop-up box summarizing the details along with a link to view the project folder on Cavatica. If the pop-up box disappears before you have a chance to click on the project link, you can [login to Cavatica](https://cavatica.sbgenomics.com){:target="_blank"} and follow the steps to [view files in Cavatica](./rna_seq_4.md#step-1-view-files-in-cavatica).

![Successful copy to Cavatica](../rna-seq-images/8_KFDRC.png "Successful copy to Cavatica")

!!! info "Query link"

    The KF Portal enables sharing of the query with the unique filter combinations including as a short URL. Login to your KF account and [click on the query link](https://p.kfdrc.org/s/6ic){:target="_blank"} to obtain the selected cohort.

    ![Sharing query](../rna-seq-images/9_KFDRC.png "Sharing query")

    [You can learn more about the different options to save/share queries in the KF Portal from our lesson](../../Common-Fund-Tools/Kids-First/Advanced-KF-Portal-Queries/KF_13_SavingQueries.md){:target="_blank"}.

In our next lesson, we will explore the newly created project folder and files on the Cavatica platform!

## Media resources

A video walkthrough of the cancer cohort selection on Kids First Portal:

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_1568tbw7&flashvars[mediaProtocol]=rtmp&amp;flashvars[streamerType]=rtmp&amp;flashvars[streamerUrl]=rtmp://www.kaltura.com:1935&amp;flashvars[rtmpFlavors]=1&amp;flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[mediaProxy.mediaPlayTo]=79&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=1_mx53rd98" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>
