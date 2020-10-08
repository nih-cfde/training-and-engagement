# Install VCFtools

VCF stands for Variant Call Format. It is a standardized text file format for representing SNP, indel, and structural variation calls. Read more about VCF files at the [GATK forum](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format).

[VCFtools](http://vcftools.sourceforge.net/man_latest.html#EXAMPLES) is a program package designed for working with VCF files. The aim of VCFtools is to provide easily accessible methods for working with complex genetic variation data in the form of VCF and BCF files. In this tutorial, you will use vcftools to convert the ".vcf" file into a format that PLINK likes.


## Step 1: Download vcftools from github

To install vcf tools in the home directory, run the following lines of code:

=== "AWS Instance"
    ```
    cd
    git clone https://github.com/vcftools/vcftools.git
    cd vcftools
    ```
The first line of code changes your directory to the home directory. The second line of code does the downloading and the last line of code takes you into the newly created directory called "vcftools".

## Step 2: Compile vcftools

Compile vcftools by running this code:

=== "AWS Instance"
    ```
    ./autogen.sh
    ./configure
    make
    sudo make install
    ```

If there are no errors at this point, you are good to go!

Look inside the "vcftools" folder to see what's in it:

=== "AWS Instance"
    ```
    ls
    ```
=== "Expected Output"

    ```
    LICENSE      README.md       build-aux    config.log     depcomp     src
    Makefile     aclocal.m4      compile      config.status  examples    stamp-h1
    Makefile.am  autogen.sh      config.h     configure      install-sh
    Makefile.in  autom4te.cache  config.h.in  configure.ac   missing
    ```
