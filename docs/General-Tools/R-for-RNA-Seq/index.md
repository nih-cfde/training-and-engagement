---
layout: page
title: Overview
---

<div class="banner"><span class="banner-text">Lesson in Development</span></div>


An Introduction to R for RNA-Seq Analysis
================================================

RNA-Sequencing (RNA-Seq) is a popular method for
determining the presence and quantity of RNA in biological samples. In
this 3 hour workshop, we will use R to explore publicly-available
RNA-Seq data from the [Gene Expression Tissue Project
(GTEx)](https://gtexportal.org/home/). Attendees will be introduced to
the R syntax, variables, functions, packages, and data structures common
to RNA-Seq projects. We will use RStudio to import, tidy, transform, and
visualize RNA-Seq count data. Attendees will learn tips and tricks for
making the processes of data wrangling and data harmonization more
manageable. This workshop will not cover cloud-based workflows for
processing RNA-seq reads or statistics and modeling because these topics
are covered in our [RNA-Seq Concepts](https://osf.io/kj5av/) and
[RNA-Seq in the
Cloud](https://training.nih-cfde.org/en/latest/Bioinformatic-Analyses/RNAseq-on-Cavatica/rna_seq_1/)
workshops. Rather, this workshop will focus on general R concepts
applied to RNA-Seq data. 

| Est. Time | Lesson name | Description |
| --- | --- | --- | 
| 25 min |[Introduction](./intro.md) | Overview of RStudio and Binder
| 30 min |[Import Data](./import.md) | Importing data with `read.csv` and `read.table`
| 30 min |[Visualize Data](./visualize.md) | Visualizing data with `ggplot2`
| 90 min |[Wrangle Data](./wrangle.md) | Tyding and transforming data
| 5 min | [Wrap-up](./wrapup.md) | Resources

### Learning Objectives

!!! info "In this workshop, you will learn how to use R and RStudio to:"

    - import and view files commonly associated with RNA-sequencing experiments
    - select variables and observations that are relevant to research questions (tidy)
    - create and rename variables (transform)
    - join data frames by common variables (harmonize)
    - visualize data using bar graphs, scatter plots, and box plots


=== "Tutorial Resources"  
    Please refer to the [RStudio cheat sheets](https://www.rstudio.com/resources/cheatsheets/) for commonly used R functions. The [R for Data Science](https://r4ds.had.co.nz/) book provides in depth descriptions and examples of many functions and concepts covered in this course.
    
    Here are the notes from workshops taught with these materials.
    
    - May 11, 2022 [workshop notes](https://hackmd.io/MsWY1O9GQXGVDl2OmV2jxg)
    - April 27, 2022 [workshop notes](https://hackmd.io/SnorsWTbTTyRenptpjrhww)
    - March 23, 2022 [workshop notes](https://hackmd.io/2ArmQGwGT0uUyL5Ehqy0hQ)

    Here is a 3 hours video of this workshop being taught at the [May 2022 Hackathon](https://nih-cfde.github.io/2022-may-hackathon/about/).

    <iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=1_61sjku7o&flashvars[localizationCode]=en&flashvars[leadWithHTML5]=true&flashvars[sideBarContainer.plugin]=true&flashvars[sideBarContainer.position]=left&flashvars[sideBarContainer.clickToClose]=true&flashvars[chapters.plugin]=true&flashvars[chapters.layout]=vertical&flashvars[chapters.thumbnailRotator]=false&flashvars[streamSelector.plugin]=true&flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&flashvars[dualScreen.plugin]=true&flashvars[Kaltura.addCrossoriginToIframe]=true&&wid=1_7ufv4tsl" width="550" height="285" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="CFDE May Hackathon - Intro to R for RNA-Seq"></iframe>. 

=== "Prerequisites"  
    Familiarity with R and RNA Sequencing is not required but would be useful. This lesson uses a standardized binder environment, which will work on Windows, Mac, and Linux operating systems, and Firefox, Safari, and Chrome web browsers.

