---
layout: page
title: Setup DESeq2 Public App
---
<style>
.highlight_txt {
  color: #9f0bde;
  background-color: #ededed;
  padding: .25em;
}
</style>

[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is a Bioconductor package used to perform DGE analysis by fitting [negative binomial model](https://www.statisticshowto.com/negative-binomial-experiment/) to the count data. It requires counts table as input along with a phenotype file describing the experimental groups.

DESeq2 performs multiple steps including

* estimating size factors which accounts for difference in library depth
* estimating gene-wise dispersions to generate accurate estimates of within-group variation
* shrinkage of dispersion estimates which reduces false positives in the DGE analysis
* perform hypothesis testing using [Wald test](https://www.statisticshowto.com/wald-test/) or [Likelihood Ratio test](https://www.statisticshowto.com/likelihood-ratio-tests/)

DESeq2 automatically removes outlier genes from analysis using [Cook's distance](https://www.statisticshowto.com/cooks-distance/) and filters genes with low counts which helps improve detection power by making the multiple testing adjustment of the p-values less severe. Refer to the [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) for detailed explanation, helpful suggestions and examples.

Cavatica offers DESeq2 as a stand alone public app which is a [Common Workflow Language (CWL)](https://www.commonwl.org) wrapper around a script with functions from the DESeq2 package. In this lesson we learn to copy, edit and setup the DESeq2 app in the project folder with cancer data files.

!!! info "Terminology"

    * Count data - represents the number of sequence reads that originated from a particular gene.
    * Dispersion - it is a measure of spread or variability in the data. DESeq2 dispersion estimates are inversely related to the mean and directly related to variance.
    * LFC - log2 fold change

## Step 1: Search & copy DESeq2 app

The first step is to obtain a copy of the DESeq2 app in the project folder.

  * Click the <span class="highlight_txt">Apps</span> tab which is currently empty and click <span class="highlight_txt">Add Apps</span> button which opens the list of Public Apps.
  * You can find the <span class="highlight_txt">DESeq2</span> app by typing any one of these keywords including "Differential Expression", "DGE", "DESeq", "RNASeq" in search bar.
  * In the <span class="highlight_txt">DESeq2</span> app box select the <span class="highlight_txt">Other versions</span> drop down box and click on the version 1.18.1. This open the app in a new tab where you can click on the `...` on the right hand corner and click <span class="highlight_txt">Copy</span>.
  * Navigate to your project Dashboard using <span class="highlight_txt">Projects</span> drop down menu and view the app under the <span class="highlight_txt">Apps</span> tab.

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_r4t7no30&flashvars[mediaProtocol]=rtmp&amp;flashvars[streamerType]=rtmp&amp;flashvars[streamerUrl]=rtmp://www.kaltura.com:1935&amp;flashvars[rtmpFlavors]=1&amp;flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[mediaProxy.mediaPlayTo]=51&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=1_cr6jqtun" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>


## Step 2: Edit DESeq2 app

The DESeq2 app has a bug with the IgnoreTxVersion parameter that can be rectified by editing the app using the tool editor.

* To do so, click on <span class="highlight_txt">DESeq2</span> in the <span class="highlight_txt">Apps</span> tab. This opens the app page.
* Click the <span class="highlight_txt">Edit</span> button on right hand upper corner which prompts a popup box with a warning message about losing update notifications for the original app. Click <span class="highlight_txt">Proceed anyway</span>.
* In the DESeq2's tool editor, find the <span class="highlight_txt">IgnoreTxVersion</span> input port and click on it.
* In the Value transform field of the port, enter the following code and click <span class="highlight_txt">Save</span>.

    ```
    {
        if ($job.inputs.ignoreTxVersion) {
          return "TRUE"
        }
       else {
          return "FALSE"
        }
    }
    ```

* Click the save icon on the top right hand corner to add a revision note.
* On the app page, the revision history is updated to read <span class="highlight_txt">Revision 1</span>.

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_8j60ve26&flashvars[mediaProtocol]=rtmp&amp;flashvars[streamerType]=rtmp&amp;flashvars[streamerUrl]=rtmp://www.kaltura.com:1935&amp;flashvars[rtmpFlavors]=1&amp;flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[mediaProxy.mediaPlayTo]=70&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=1_x2wd8enj" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>


## Step 3: Obtain reference gene annotation

A reference gene annotation file in GTF format is required by DESeq2 app to summarize the transcript level abundances contained in the [Kallisto](http://pachterlab.github.io/kallisto//releases/2017/03/20/v0.43.1) files for gene-level analysis. Internally, [tximport](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html) another Bioconductor package, is utilized to obtain the gene level summary.  

* Edit the metadata columns to show Reference genome under <span class="highlight_txt">Files</span> tab. All the files used the GRCh38 (hg38) homo sapiens genome assembly released by Genome Reference Consortium.
* Click on <span class="highlight_txt">Data</span> drop down menu and click on <span class="highlight_txt">Public Reference Files</span>.
* This takes you to a new page for <span class="highlight_txt">Public Files</span>. Click on <span class="highlight_txt">Type: All</span> button to bring a drop down list and select **GTF**.
* From the results, select the ENSEMBL Release 84 version of the Human gene annotation in GTF format. Click on <span class="highlight_txt">Copy</span> and select the project folder with the cancer files. Select <span class="highlight_txt">Copy</span> in the popup window.
* A notification menu will highlight the successful copy of the file and clicking on the project folder name will take you to the <span class="highlight_txt">Files</span> tab in folder.
* Check for the reference file using the <span class="highlight_txt">Type: All</span> button and select **GTF**.

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_pxn0zxoc&flashvars[mediaProtocol]=rtmp&amp;flashvars[streamerType]=rtmp&amp;flashvars[streamerUrl]=rtmp://www.kaltura.com:1935&amp;flashvars[rtmpFlavors]=1&amp;flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[mediaProxy.mediaPlayTo]=48&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=1_l7cbi5kv" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>
