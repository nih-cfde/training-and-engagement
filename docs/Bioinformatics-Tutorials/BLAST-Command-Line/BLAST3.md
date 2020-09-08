---
layout: page
title: Install BLAST
---

## Set-up an Instance

Boot a t2.micro instance on AWS and connect your shell prompt.

!!! Tip 
	
	If you need assistance setting up an instance check out the tutorial [Intro to AWS](<../../General-Tutorials/Introduction_to_Amazon_Web_Services/introtoaws1.md>)!


# Install BLAST software using Command Line
Now, install some software. We will need NCBI BLAST for the below tutorial


Copy and paste the following commands
```
   sudo apt-get update && sudo apt-get -y install python ncbi-blast+
```

sudo apt-get update command is used to download package information from all configured sources

So when you run update command, it downloads the package information from the internet. This updates the software list and installs the Python programming
language and NCBI BLAST+.
