---
layout: page
title: Install PLINK
---

Install PLINK
===============

## What is PLINK?
[PLINK](http://zzz.bwh.harvard.edu/plink/index.shtml) is a free, open-source whole-genome association analysis toolset. The software is designed flexibly to perform a wide range of basic, large-scale genetic analyses including GWAS.

## Installation

To install PLINK, type these commands into your [Ubuntu AWS terminal](./Accessing_aws.md):

```
cd /usr/local/bin/
sudo wget https://zzz.bwh.harvard.edu/plink/dist/plink-1.07-x86_64.zip

```

!!! Error
    If you see a `wget` error message when you run this code, you will have to install `wget` by running these lines of code:


    `/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`

    `brew install wget`


Next, unzip and cd into the folder:

```
sudo unzip -o plink-1.07-x86_64.zip
sudo rm plink-1.07-x86_64.zip
cd plink-1.07-x86_64
```

!!! Error
    If you see a error like this:

    `sudo: unzip: command not found`

    You will need to install unzip:

    `sudo apt-get install unzip`

    Then re-run the unzip code block.

After cd-ing into the plink-1.07-x86_64 folder, add plink to the `.bashrc`. By doing so, you are telling the computer that PLINK can be accessed from any directory. Type:

```
echo export PATH=$PATH:$(pwd) >> ~/.bashrc
source ~/.bashrc
```

PLINK has now been successfully installed!
