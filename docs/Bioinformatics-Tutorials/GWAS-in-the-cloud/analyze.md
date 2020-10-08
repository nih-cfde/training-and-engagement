---
layout: page
title: Analyze
---

Summary Statistics and Association Analysis
===========================================

All of these steps must be performed on your Ubuntu AWS terminal window.

## Step 1: Convert VCF into PLINK readable format

Remember that VCF files are variant calling format files that have a very [specific structure](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format). PLINK does not take vcf files as inputs. So you must convert the ".vcf" into PLINK readable format: ped and map.

[PED and MAP files](http://zzz.bwh.harvard.edu/plink/data.shtml) are plain text files; ped files contain genotype information (one individual per row) and map files contain information on the name and position of the markers in the PED file.

First change directories to "GWAS":

=== "AWS Instance"
    ```
    cd /home/ubuntu/GWAS
    ```

Then make the map and ped files:

=== "AWS Instance"
    ```
    vcftools --vcf pruned_coatColor_maf_geno.vcf --plink --out coatColor
    ```

the --plink options outputs the genotype data in PLINK ped format. Two files are generated, with suffixes ".ped" and ".map"

!!! Error
    If you get a vcftools install error, follow the directions in the error message to install vcftools. Visit the [vcftools](./vcftools_install.md) page of this tutorial for detailed installation instructions.

## Step 2: Create list of minor alleles

!!! note "Population Genetics Terms"

    **For a given locus**

    Major allele: is the most common allele in the population

    Minor allele: is the least common allele in the population

    Risk allele: in the context of a disease, is the allele associated with the disease.
    For most Mendelian diseases, and a few (multi-gene) complex diseases, the risk allele is the minor allele. However, in some cases, the risk allele can be the major allele.

    For the purposes of this tutorial, we will set the minor allele at each SNP locus to be the risk (or reference) allele. This makes visualization and interpretation of results easier.


In order to specify the minor allele as the reference allele for PLINK (A1), you must create a list of these alleles. We're calling the list "minor_alleles". To do so, run:

=== "AWS Instance"
    ```
    cat pruned_coatColor_maf_geno.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{{if($3==".")$3=$1":"$2;}print $3,$5;}'  > minor_alleles
    ```

A detailed explanation of [`awk`](https://www.grymoire.com/Unix/Awk.html) is beyond the scope of this tutorial. However, the gist of this code is that it grabs the vcf file, extracts the third (SNP position info) and fifth (minor allele info) columns, and outputs them into a file called "minor_alleles".

!!! Note
    In the vcf file, the REF column contains major alleles and the ALT column contains minor alleles. For this tutorial, we are grabbing the ALT column from the vcf file and using it to set the minor alleles as the reference alleles in PLINK.


## Step 3: Quality Control

Quality control (QC) is an important step in GWAS and must be done per individual and per marker.

**Per-individual** QC of GWA data consists of identifying individuals with: 1) missing genotype or heterozygosity rate. 2) replicated samples or closely related individuals. 3) identification of individuals of divergent ancestry.

**Per-marker** QC of GWA data consists of identifying SNPs with: 1) lots of missing genotypes. 2) significant deviation from Hardy-Weinberg equilibrium (HWE). 3) largely different missing genotype rates between cases and controls. 4) very low minor allele frequencies

Read more about quality control in this journal article by [Anderson et al. 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3025522/)

### Missing rates
In this tutorial, we will generate some simple summary statistics on rates of missing data in the file, using the [--missing option](http://www.cog-genomics.org/plink/1.9/basic_stats#missing):

=== "AWS Instance"
    ```
    plink --file coatColor --make-pheno coatColor.pheno "yellow" --missing --out miss_stat --noweb --dog --reference-allele minor_alleles --allow-no-sex --adjust
    ```
=== "Expected Output"

    ```
    @----------------------------------------------------------@
    |        PLINK!       |     v1.07      |   10/Aug/2009     |
    |----------------------------------------------------------|
    |  (C) 2009 Shaun Purcell, GNU General Public License, v2  |
    |----------------------------------------------------------|
    |  For documentation, citation & bug-report instructions:  |
    |        http://pngu.mgh.harvard.edu/purcell/plink/        |
    @----------------------------------------------------------@

    Skipping web check... [ --noweb ]
    Writing this text to log file [ miss_stat.log ]
    Analysis started: Wed Sep 16 20:46:53 2020

    Options in effect:
	   --file coatColor
	   --make-pheno coatColor.pheno yellow
	   --missing
	   --out miss_stat
	   --noweb
	   --dog
	   --reference-allele minor_alleles
	   --allow-no-sex
	   --adjust

     476840 (of 476840) markers to be included from [ coatColor.map ]
     Warning, found 53 individuals with ambiguous sex codes
     Writing list of these individuals to [ miss_stat.nosex ]
     53 individuals read from [ coatColor.ped ]
     0 individuals with nonmissing phenotypes
     Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
     Missing phenotype value is also -9
     0 cases, 0 controls and 53 missing
     0 males, 0 females, and 53 of unspecified sex
     Constructing a binary phenotype from [ coatColor.pheno ]
     Test value is [ yellow ] and missing value is [ -9 ]
     53 of 53 individuals assigned to 2 cluster(s)
     Set 24 cases and 29 controls, 0 missing, 0 not found
     Reading SNPs to set reference allele [ minor_alleles ]
     Set reference alleles for 476840 SNPs, 424930 different from minor allele
     Before frequency and genotyping pruning, there are 476840 SNPs
     53 founders and 0 non-founders found
     Writing individual missingness information to [ miss_stat.imiss ]
     Writing locus missingness information to [ miss_stat.lmiss ]
     Total genotyping rate in remaining individuals is 0.977191
     0 SNPs failed missingness test ( GENO > 1 )
     0 SNPs failed frequency test ( MAF < 0 )
     After frequency and genotyping pruning, there are 476840 SNPs
     After filtering, 24 cases, 29 controls and 0 missing
     After filtering, 0 males, 0 females, and 53 of unspecified sex

     Analysis finished: Wed Sep 16 20:47:10 2020
    ```

!!! Note
        **What are all these PLINK tags?**

        `--file`: tells it the name of PLINK readable files

        `--missing`: produces sample-based and variant-based missing data reports using default filters

        `--out`: name of the output file

        `--dog`: tells PLINK to look at the dog genome
        The default reference genome option is human. Other available options are: `--mouse`, `--horse`, `--cow` and `--sheep`

        `--make-pheno`: tells PLINK to look at the coatColor.pheno file for phenotype information and sets the alternative phenotype to "yellow"

        `--reference-allele`: sets the A1 or minor allele using the file minor_alleles

        `--allow-no-sex`: since our dataset does NOT have a "sex" field, this option allows plink to ignore the missing sex field

        `--noweb`: each time PLINK runs, it checks for an update. On a slow network this sometimes causes delays and the --noweb option disables this


There's a lot of information in the output, but the relevant bits:

```
476840 (of 476840) markers to be included from [coatColor.map]
```
Indicates that all markers can be included.

```
53 individuals with ambiguous sex codes
```
There is no column for sex in our dataset. That's fine, you told PLINK to ignore sex.

```
Test value is [yellow] and missing value is [-9]
53 of 53 individuals assigned to 2 cluster(s)
Set 24 cases and 29 controls, 0 missing, 0 not found
```
This means 53 individuals have been assigned to two clusters: test (yellow coat color) and control (dark coat color). There are 24 yellow coat color and 29 dark coat color individuals, and no individuals have missing phenotype data.

```
Total genotyping rate in remaining individuals is 0.977
```
About 2% of genotypes are missing after thresholding.

```
0 SNPs failed missingness test (GENO>1)
0 SNPs failed frequency test (MAF<0)
```
Here, GENO>1 means exclude an individual if all of its genotypes are missing. Obviously, this is a pretty lenient parameter. Similarly, MAF (minor allele frequency)<0 means exclude all minor alleles that have a frequency lower than 0. You may wish to change these thresholds based on your research question by explicitly specifying `--mind` or `--geno` or `--maf`.

The per individual and per SNP rates are then output to the files "miss_stat.imiss" and "miss_stat.lmiss", respectively. If you had not specified an `--out` option, the root output filename would have defaulted to "plink".

Look at the per SNP rates by running:

=== "AWS Instance"
    ```
    less miss_stat.lmiss
    ```
=== "Expected Output"
    ```
    CHR                                    SNP   N_MISS   N_GENO   F_MISS
    1                          BICF2P1489653        1       53  0.01887
    1                             chr1:11368        0       53        0
    1                        BICF2G630707787        0       53        0
    1                             chr1:22137        0       53        0
    1                             chr1:22143        0       53        0
    1                             chr1:23623        0       53        0
    1                             chr1:23651        2       53  0.03774
    1                             chr1:23653        2       53  0.03774
    1                         BICF2S23441188        1       53  0.01887
    1                             chr1:30122        2       53  0.03774
    1                             chr1:30173        2       53  0.03774
    1                             chr1:30370        1       53  0.01887
    1                             chr1:30664        0       53        0
    1                             chr1:30723        0       53        0
    1                             chr1:31554        1       53  0.01887
    1                             chr1:32439        1       53  0.01887
    1                             chr1:33435        1       53  0.01887
    1                             chr1:34323        0       53        0
    1                           BICF2P476457        1       53  0.01887
    1                        BICF2G630707798        0       53        0
    1                             chr1:35781        1       53  0.01887
    1                              rs8471230        0       53        0
    ```


For each SNP, you see the number of missing individuals (N_MISS) and the proportion of individuals missing (F_MISS).
For examples, the SNP `BICF2P1489653` is missing in 1 out of 53 individuals, giving it a missing frequency of 0.01886792452 (i.e. 1/53). Lower proportions are better!

You can quit this mode and return to the terminal by typing ++q++.

Similarly, look at the per individual rates in the "miss_stat.imiss" by typing

=== "AWS Instance"
    ```
    less miss_stat.imiss
    ```
    You can quit this mode and return to the terminal by typing ++q++.

=== "Expected Output"
    ```
      FID       IID MISS_PHENO   N_MISS   N_GENO   F_MISS
      dark_13   dark_13          N     4994   476840  0.01047
      dark_23   dark_23          N     4478   476840 0.009391
      dark_21   dark_21          N     4739   476840 0.009938
      yellow_5  yellow_5         N    15094   476840  0.03165
      yellow_6  yellow_6         N    13889   476840  0.02913
      dark_1    dark_1           N     5703   476840  0.01196
      dark_7    dark_7           N     4895   476840  0.01027
      dark_9    dark_9           N    23771   476840  0.04985
      dark_2    dark_2           N    33736   476840  0.07075
      dark_10   dark_10          N    20475   476840  0.04294
      dark_8    dark_8           N     5331   476840  0.01118
      yellow_24 yellow_24        N     3932   476840 0.008246
      yellow_22 yellow_22        N     5199   476840   0.0109
      dark_17   dark_17          N      448   476840 0.0009395
      dark_18   dark_18          N     1622   476840 0.003402
      dark_4    dark_4           N     5691   476840  0.01193
      dark_6    dark_6           N     5419   476840  0.01136
      dark_5    dark_5           N    18523   476840  0.03885
      yellow_16 yellow_16        N     7534   476840   0.0158
      yellow_14 yellow_14        N     1640   476840 0.003439
      yellow_18 yellow_18        N     1657   476840 0.003475
      yellow_15 yellow_15        N     4159   476840 0.008722
    ```

The final column is the genotyping rate for that individual. Looking at the first row, the individual dark_13 has 4994 missing SNPs out of 476840, producing a missing genotype rate of 0.01047.

In this tutorial, we are *not* excluding any SNPs or individuals from downstream association analyses. However, if the missing genotype rate per SNP or individual is high, PLINK has tags to exclude those genotypes or individuals based on [user-specified criteria](http://www.cog-genomics.org/plink/1.9/filter).

## Step 4: Convert to PLINK binary format

Next, convert the output file (coatColor) to PLINK binary format (fam,bed,bim) for downstream analysis:

=== "AWS Instance"
    ```
    plink --file coatColor --allow-no-sex --dog --make-bed --noweb --out coatColor.binary

    ```

[`--make-bed`](http://www.cog-genomics.org/plink/1.9/data) creates a new PLINK binary file set, after applying sample/variant filters and other operations.



## Step 5: Run a simple association analysis

Learn more about [association tests here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002822#s7)

=== "AWS Instance"

    ```
    plink --bfile coatColor.binary --make-pheno coatColor.pheno "yellow" --assoc --reference-allele minor_alleles --allow-no-sex --adjust --dog --noweb --out coatColor
    ```

!!! Note
    **What are these new PLINK tags?**

    `--bfile`: takes .binary file as input.

    `--assoc`: performs a standard case/control association analysis which is a chi-square test of allele frequency.

    `--adjust`: enables correction for multiple analysis and automatically calculates the genomic inflation factor
