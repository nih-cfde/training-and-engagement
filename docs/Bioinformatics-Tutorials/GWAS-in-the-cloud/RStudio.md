---
layout: page
title: Install R and RStudio
---

Install R and RStudio
=====================


## What is R and RStudio?

R is a language and environment for statistical computing and graphics. [R studio](https://rstudio.com/) is an integrated development environment for R. Visit [this website](https://rstudio.com/products/rstudio/download-server/debian-ubuntu/) for more details on downloading.

You will use R to visualize loci that are strongly associated with coat color.

## Installation

To download R and RStudio, type:

```
sudo apt-get update && sudo apt-get install -y gdebi-core r-base r-base-dev
wget https://download2.rstudio.org/server/xenial/amd64/rstudio-server-1.3.959-amd64.deb
sudo gdebi rstudio-server-1.3.959-amd64.deb
```
!!! Error

    If you get an error message about gdebi, you will need to install it:

    `sudo apt-get install gdebi-core`

    Then re-run
    `sudo gdebi rstudio-server-1.3.959-amd64.deb`
