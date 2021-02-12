---
layout: page
title: Analysis using Data Cruncher
---

So far we have explored running DGE analysis using a public app based on DESeq2. In the second approach, we will setup an interactive analysis on an instance running RStudio computational environment. We will run DGE workflow using an analysis script, generate reports and plots.

!!! note "DGE Tools"

    While there are other established tools to perform DGE analysis including [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and [Limma-Voom](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29), we will be using DESeq2 in our script to allow you to compare the output between the two approaches.


## Step 1: Starting Data Cruncher

* Click the <span class="highlight_txt">Interactive Analysis</span> tab located on the right hand corner below your account settings menu.
* Select <span class="highlight_txt">Open</span> in the Data Cruncher panel.
* Click on the <span class="highlight_txt">Create your first analysis</span> which appears the first time your are setting up.
* In the popup box, select <span class="highlight_txt">RStudio</span> for Environment. Provide an analysis name in the box. Here, **Cancer_DGE** was used to title the analysis. Click <span class="highlight_txt">Next</span> when done.
* In the `Compute requirements` tab, we will use the default instance type (c5.2xlarge, $0.49/hr). We increase the `Suspend time`, which is the period of inactivity after which the instance is stopped automatically and the analysis is saved, from 30 to 60 minutes.
* Click <span class="highlight_txt">Start the analysis</span>. This prompts the initializing of the analysis which involves set up of the instance and  preparing analysis environment.

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_4hc4cd8i&flashvars[mediaProtocol]=rtmp&amp;flashvars[streamerType]=rtmp&amp;flashvars[streamerUrl]=rtmp://www.kaltura.com:1935&amp;flashvars[rtmpFlavors]=1&amp;flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[mediaProxy.mediaPlayTo]=44&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=1_6ih50uzu" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>

!!! info "Instance Types"

    You can find details on all available US instances from Amazon Web Services (AWS) on [Cavatica's Platform Documentation](https://docs.sevenbridges.com/docs/list-of-available-amazon-web-services-instances).

## Step 2: Navigating analysis editor and load script

After the instance is initialized, you will be automatically directed to the analysis editor which in this case is the RStudio interface.

!!! info "RStudio IDE"

    [Read more about the different panes and options of the RStudio interface](https://georgejmount.com/tourofrstudio/), the integrated development environment (IDE) for the R programing language.

### Directory structure

The editor is associated with a directory structure to help you navigate the working space. You can access it via the Files/Packages/Plots/Help/Viewer pane on the bottom right hand corner of RStudio.

```
/sbgenomics
|–– output-files
|–– project-files
|–– projects
|–– workspace
```

!!! important "Important"

    The `project-files` directory which contains all the input files is a read only file system while you have read-write permissions for the `workspace` and `output-files` directories.

* **workspace** is the default working directory for the analysis. You can use the RStudio <span class="highlight_txt">Upload</span> option to get files from your local computer to the workspace.
* **output-files** can be used as the directory to save all the outputs from your analysis. If not specified, the files are saved to workspace
* **project-files** is the directory containing all the input files from the current project. Since it is a read only file system, no changes can be made to these files via the editor interface.

You can access the scripts for this lesson in [Step 3](#R-code). To upload the script file, <a name="Upload-file"></a>

 * Click on the <span class="highlight_txt">Upload</span> option in the Files/Packages/Plots/Help/Viewer pane.
 * Click <span class="highlight_txt">Choose File</span> to select the file from your local computer.
 * Once uploaded, click on the script file name to open it in the script editor pane (top left hand corner).

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_uennaoay&flashvars[mediaProtocol]=rtmp&amp;flashvars[streamerType]=rtmp&amp;flashvars[streamerUrl]=rtmp://www.kaltura.com:1935&amp;flashvars[rtmpFlavors]=1&amp;flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[mediaProxy.mediaPlayTo]=33&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=1_77fny8c6" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>

### Session outputs

The generated output and environment files from an active session are saved when the analysis is stopped by clicking <span class="highlight_txt">Stop</span> located on the right hand top corner. You can access the session files via the <span class="highlight_txt">Files</span> tab in your project folder.

The Data Cruncher comes with set of libraries that are pre-installed. These vary depending on the environment you chose during setup. We chose the default environment for RStudio `SB Bioinformatics - R 4.0` which is loaded with [set of CRAN and Bioconductor libraries](https://docs.sevenbridges.com/docs/about-libraries-in-a-data-cruncher-analysis).

!!! important "Installing additional libraries"

    Although the output files, the environment and history of the session are saved upon stopping the analysis editor, any installed libraries are only good for the session and must be re-installed for every restart of the instance.


## Step 3: Run analysis script <a name="R-code"></a>

 We provide you with the option to download two versions of the analysis script based on choice of execution in RStudio: </br>
 (a) [version to execute automatically using `Source`](./rna-seq-supporting-docs/Cancer_DGE_Analysis_Automate.R) </br>
 (b) [version to execute the code in chunks using the `Run` option](./rna-seq-supporting-docs/Cancer_DGE_Analysis.R).

You can review the steps to [**upload**](#Upload-file) the R script to the **workspace** directory.
The version of the script that will be run manually contains some additional packages and lines of code to allow for interactive exploration of the data prior to analysis. The DGE analysis and all the generated output are otherwise identical between the two versions.   

!!! important "Phenotype File Name"

      For the scripts to run error free ensure that the name of the phenotype CSV file is [**phenotype_filtered.csv**](./rna_seq_6.md#step-3-upload-phenotype-file-to-cavatica). Update the scripts to a different name if applicable before execution.

### Step 3a: Execute using `Source`

To get started click on the down arrow next to <span class="highlight_txt">Source</span> and click <span class="highlight_txt">Source with Echo</span>. This will print the comments as the code is executed.

![Source with Echo](../rna-seq-images/rna-seq-8-1.png "Source with Echo")

Once completed, you will get a popup window asking to try to open the html report. Click </br><span class="highlight_txt">Try Again</span> to open a new tab for the report.

![Popup window](../rna-seq-images/rna-seq-8-2.png "Popup window")

Alternatively, you can click <span class="highlight_txt">Cancel</span> in the popup window and subsequently click <span class="highlight_txt">Stop</span> to view the files in your project folder.

![Stop analysis](../rna-seq-images/rna-seq-8-3.png "Stop analysis")

For costs and time comparison between the two approaches, we use the automated version with the option to view the output files in the project folder which took 25 minutes to run and cost $0.2.


### Step 3b: Execute using `Run`

You can also execute the code by selecting chunks or line by line and clicking the <span class="highlight_txt">Run</span> option or using ++ctrl+enter++ keys. This allows you greater flexibility to explore and understand the data. The first step is installing the packages necessary for DGE analysis and this takes approximately 17 minutes.

![Use Run](../rna-seq-images/rna-seq-8-4.png "Use Run")

This version include the Bioconductor package [`pcaExplorer`](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2879-1), which provides interactive visualization of RNA-Seq datasets based on [Principal Components Analysis](https://www.youtube.com/watch?v=_UVHneBUBW0&list=PLblh5JKOoLUJo2Q6xK4tZElbIvAACEykp&index=22). Watch the video below to learn how to use `pcaExplorer` for the filtered cancer dataset.

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_1z57clip&flashvars[mediaProtocol]=rtmp&amp;flashvars[streamerType]=rtmp&amp;flashvars[streamerUrl]=rtmp://www.kaltura.com:1935&amp;flashvars[rtmpFlavors]=1&amp;flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[mediaProxy.mediaPlayTo]=104&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=1_tywvb439" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>



After having run the entire script, click <span class="highlight_txt">Stop</span> to view the files in your project folder.

!!! bug "Login Timeout"

    It is possible to be logged out of Cavatica despite having an active RStudio session. You will be unable to stop the analysis from within the editor using <span class="highlight_txt">Stop</span> if that occurs.

     * Login to Cavatica in a new tab or window.
     * Navigate to the data cruncher session via either the <span class="highlight_txt">Interactive Analysis</span> tab or using the `ANALYSES` pane in your project home page.

     ![Data Cruncher quick access](../rna-seq-images/rna-seq-8-5.png "Data Cruncher quick access"){: width=70%}

    * Click <span class="highlight_txt">Stop</span> on the session page.

## Step 4: View output files

All the session files and the generated outputs are saved after the analysis is stopped and are accessible to you on the session page.

![Output Files](../rna-seq-images/rna-seq-8-6.png "Output Files")

The tag for the session changes from {==RUNNING==} to {==SAVED==}. Similar to the DESeq2 app, four output files are generated:

* **Cancer_DESeq2_DGE_results.csv** contains the ordered table of gene level statistics generated using the `results()` function in DESeq2 package.
* **Cancer_DESeq2_normalized_counts.txt** contains counts normalized using the estimated sample-specific normalization factors.
* **DESeq2-Report** folder which contains the HTML report generated using [regionReport](https://f1000research.com/articles/4-105/v2). The report contains all the visualizations along with the associated code from the [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#how-do-i-use-vst-or-rlog-data-for-differential-testing).
* **Cancer_DGE_{Date}.env.RData** is the R workspace image that includes all the objects and variables generated from the code. The `.RData` listed under Workspace is saved by default by the Data Cruncher.

!!! info "Output Differences"

        Although the DGE results are the same between the two analysis approaches, there are some differences between the two `html` reports since they are not the exact code. The MA plots generated using Data Cruncher use blue to signify significant genes and the counts plot use points instead of bars.

All the files are clickable for preview on Cavatica. You can either download individual files by clicking on the file name or follow the steps to [tag and download the files listed in the Analysis with DESeq2 Public App](./rna_seq_7.md#step-4-tag-download-analysis-outputs) lesson.

## Conclusion

This concludes the RNAseq on Cavatica tutorial. We hope that you found the tutorial helpful and will continue to use cloud computing for your analysis!

!!! note " Key Points"

    * [Kids First Portal](https://kidsfirstdrc.org) is the go to resource for pediatric cancer & structural birth defects datasets.
    * Data to analysis in one step using [Cavatica](cavatica.sbgenomics.com), the cloud based analysis platform integrated into Kids First Portal.
    * You can filter, view, and download data from Cavatica.
    * Upload data to Cavatica from multiple sources including your local machine.
    * You can search, copy and modify a public app on Cavatica.
    * Setup and successful run of the DESeq2 app by choosing appropriate inputs.
    * Setup a virtual computational environment running RStudio and analyze by executing code from a script.
