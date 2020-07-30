---
layout: page
title: Manhattan Plots
---

Create a Manhattan Plot
===========================================

## Install qqman package

```
sudo Rscript -e "install.packages('qqman',  contriburl=contrib.url('http://cran.r-project.org/'))"
```

## Identify statistical cutoffs

This code finds the equivalent of 0.05 and 0.01 p value in the negative-log-transformed p values file. We will use these cutoffs to draw horizontal lines in the Manhattan plot for visualization of haplotypes that cross the 0.05 and 0.01 statistical threshold (i.e. have a statistically significant association with yellow coat color)

```
unad_cutoff_sug=$(tail -n+2 coatColor.assoc.adjusted | awk '$10>=0.05' | head -n1 | awk '{print $3}')
unad_cutoff_conf=$(tail -n+2 coatColor.assoc.adjusted | awk '$10>=0.01' | head -n1 | awk '{print $3}')
```

## Install GhostScript by running:

```
sudo apt-get update -y
sudo apt-get install -y ghostscript-x
```

## Run the plotting function

```
Rscript -e 'args=(commandArgs(TRUE));library(qqman);'\
'data=read.table("coatColor.assoc", header=TRUE); data=data[!is.na(data$P),];'\
'bitmap("coatColor_man.bmp", width=20, height=10);'\
'manhattan(data, p = "P", col = c("blue4", "orange3"),'\
'suggestiveline = 12,'\
'genomewideline = 15,'\
'chrlabs = c(1:38, "X"), annotateTop=TRUE, cex = 1.2);'\
'graphics.off();' $unad_cutoff_sug $unad_cutoff_conf
```

The output should look like this:

![](images/Final_Output.png)

Do an `ls` to check if you've created the .bmp file.... and it worked!

![](images/BMP.png)

## Visualization

You can visualize the Manhattan plot by downloading the coatColor_man.bmp file to our local computers. To do this, open a new terminal window like before (i.e. by selecting the terminal window and typing `cmd+N`).

Now run the following code to do the actual copying:

```
scp -i ~/Desktop/amazon.pem -r ubuntu@ec2-??-???-???-??.us-east-2.compute.amazonaws.com:/home/ubuntu/GWAS/ ~/Desktop/
```
!!! Note
    The -r flag used in the copy command, copies the entire "GWAS" folder. To copy just the .bmp file, drop the -r flag and include the filename in the aws path.


The file has been copied to your Desktop! Here's what it should look like when opened using the "Preview" application in Mac:

![](images/coatColor_man.png)

The X- axis represents haplotypes from each region of the genome that was tested, organized by chromosome. Each colored block represents a chromosome and is made of thousands of dots that represent haplotypes. The Y axis is a p value (probability that the association was observed by chance) and is negative log transformed.

In our graph, halptypes in four parts of the genome (chromosome 2, 5, 28 and X) are found to be associated with an increased occurrence of the yellow coat color phenotype.

The top associated mutation is a nonsense SNP in MC1R (c.916C>T) known to control pigment production.

[Watch this video for a detailed explanation of GWAS Manhatten plots](https://www.google.com/search?q=how+to+read+a+manhattan+plot&oq=how+to+read+a+manhattan+plot&aqs=chrome..69i57.7911j0j4&sourceid=chrome&ie=UTF-8#kpvalbx=_tXIPX9mmFsmT0PEP64-OkAk26)
