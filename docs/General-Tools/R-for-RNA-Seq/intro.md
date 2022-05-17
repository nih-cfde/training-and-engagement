---
title: Introduction
---

The book [“R for Data Science”](https://r4ds.had.co.nz/index.html) provides an excellent framework for using data science to turn raw data
into understanding, insight, and knowledge. We will use this framework
as an outline for this workshop.

**R** is a statistical computing and data visualization programming
language. **RStudio** is an integrated development environment, or IDE,
for R programming. R and RStudio work on Mac, Linux, and Windows
operating systems. The RStudio layout displays lots of useful
information and allows us to run R code from a script and view and save
outputs all from one interface.

When you start RStudio, you’ll see two key regions in the interface: the
console and the output. When working in R, you can type directly into
the console, or you can type into a script. Saving commands in a script
will make it easier to reproduce. You will learn more as we go along!
  
![](./images/RStudio.png)

For today’s lesson, we will focus on data from the [Gene-Tissue
                                                    Expression (GTEx) Project](https://commonfund.nih.gov/gtex). GTEx is an
ongoing effort to build a comprehensive public resource to study
tissue-specific gene expression and regulation. Samples were collected
from 54 non-diseased tissue sites across nearly 1000 individuals,
primarily for molecular assays including WGS, WES, and RNA-Seq. 


### Getting Started

1.  Click the
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nih-cfde/training-rstudio-binder/data?urlpath=rstudio)
button to generate a computing environment for this workshop.
2.  Navigate to the GTEx folder.
3.  Click `GTEx.Rproj` and click “Yes” to open up an Rproject. This will
set the working directory to `~/GTEx/`.
4.  If you open the `r4rnaseq-workshop.R` file which contains all the
commands for today’s workshop, you can click through this and all
the commands should run successfully.
5.  If you open a new R Script by clicking **File > New File > R
Script**, you can code along by typing out all the commands for
today’s lesson as I type them.

Click “Run” to send commands from a script to the console or click
command enter.

### R is a calculator

You can perform simple and advanced calculations in R.

``` r
2 + 2 * 100
```

    ## [1] 202

``` r
log10(0.05)
```

    ## [1] -1.30103

You can save variable and recall them later.

``` r
pval <- 0.05
pval
```

    ## [1] 0.05

``` r
-log10(pval)
```

    ## [1] 1.30103

You can save really long lists of things with a short, descriptive names
that are easy to recall later.

``` r
favorite_genes <- c("BRCA1", "JUN",  "GNRH1", "TH", "AR")
favorite_genes
```

    ## [1] "BRCA1" "JUN"   "GNRH1" "TH"    "AR"

### Loading R packages

Many of the functions we will use are pre-installed. The
[Tidyverse](https://www.tidyverse.org/) is a collection of R packages
that include functions, data, and documentation that provide more tools
and capabilities when using R. You can install the popular data
visualization package `ggplot2` with the command
`install.packages("ggplot2")`). It is a good idea to “comment out” this
line of code by adding a `#` at the beginning so that you don’t
re-install the package every time you run the script. For this workshop,
the packages listed in the `.binder/environment.yml` file were
pre-installed with Conda.

``` r
#install.packages("ggplot2")
```

After installing packages, we need to load the functions and tools we
want to use from the package with the `library()` command. Let’s load
the `ggplot2` package.

``` r
library(ggplot2)
```

Now you have successfully loaded the necessary R packages. Let's complete and exercise

### Exercise

=== "Challenge"

    We will also use functions from the packages `tidyr` and `dplyr` to tidy and transform data. What command would you run to load these packages?
  
=== "Answer"

    ```r
    library(tidyr)  
    library(dplyr)
    ```

You can also navigate to the “Packages” tab in the bottom right pane of
RStudio to view a list of available packages. Packages with a checked
box next to them have been successfully loaded. You can click a box to
load installed packages. Clicking the “Help” Tab will provide a quick
description of the package and its functions.


### Key functions

| Function             | Description                                  |
|----------------------|----------------------------------------------|
| `<-`                 | The assignment variable                      |
| `log10()`            | A built-in function for a log transformation |
| `install.packages()` | An R function to install packages            |
| `library()`          | The command used to load installed packages  |