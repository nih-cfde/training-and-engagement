---
layout: page
title: Download the data
---

Download the data
=================

Make a directory for the GWAS analysis called `GWAS` and cd into that directory.

```
mkdir ~/GWAS && cd ~/GWAS
```
!!! Note
    The coat color data lives in a website called cyverse. It is not easy to make AWS talk to cyverse, so download the data onto your LOCAL computer and then upload it to AWS.

## Download onto your laptop

To download the data onto you laptop, run the following code on a **new terminal window** on your local machine. Start a new terminal by selecting the open terminal window and typing `cmd+N`. Then run:

```
wget https://de.cyverse.org/dl/d/E0A502CC-F806-4857-9C3A-BAEAA0CCC694/pruned_coatColor_maf_geno.vcf.gz
wget https://de.cyverse.org/dl/d/3B5C1853-C092-488C-8C2F-CE6E8526E96B/coatColor.pheno
```
this download will take a few seconds.

## Upload onto AWS
From the same terminal window, upload the two files onto the AWS computer. Be sure to change the path to point to your two newly downloaded files. Also change the ec2 instance name like you did before.

```
scp -i ~/Desktop/amazon.pem path/to/file/pruned_coatColor_maf_geno.vcf.gz ubuntu@ec2-???-???-???-???.compute-1.amazonaws.com:~/GWAS/
gunzip pruned_coatColor_maf_geno.vcf.gz #unzip the vcf
scp -i ~/Desktop/amazon.pem path/to/file/coatColor.pheno ubuntu@ec2-???-???-???-???.compute-1.amazonaws.com:~/GWAS/
```
