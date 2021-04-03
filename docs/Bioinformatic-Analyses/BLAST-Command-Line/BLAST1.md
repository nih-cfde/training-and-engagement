---
layout: page
title: BLAST Overview
hide:
  - toc
---


Running Command-Line BLAST
=============================


BLAST is the **B**asic **L**ocal **A**lignment **S**earch large sequence databases; It starts by finding small matches between the two sequences and extending those matches.  For in-depth information on how BLAST works and the different BLAST functionality, check out the [resources page](https://blast.ncbi.nlm.nih.gov/Blast.cgi).

BLAST can be helpful for identifying the source of a sequence or finding a similar sequence in another organism.  In this lesson, we will use BLAST to find zebrafish proteins that are similar to a small set of mouse proteins.

Why use the command line?
BLAST has a very nice graphical interface for searching sequences in NCBI's database.
However, running BLAST through the commmand line has many benefits:

- It's much easier to run many BLAST queries using the command line than the GUI

- Running BLAST with the command line is reproducible and can be documented in a script

- The results can be saved in a machine-readable format that can be analyzed later

- You can create your own databases to search rather than using NCBI's pre-built databases

- It allows the queries to be automated

- It allows you to use a remote computer to run the BLAST queries

Est. time | Lesson name | Description
--- | --- | ---
30 mins | [Install BLAST](../BLAST-Command-Line/BLAST3.md) | Set up local BLAST installation on your computer
15 mins | [How to Run BLAST+](../BLAST-Command-Line/BLAST4.md) | Run BLAST analysis with AWS

!!! note "Learning Objectives"

    - Gain hands-on exposure to the linux command line

    - Understand how data is turned into results by programs run at the command line

=== "Prerequisites"

    - Some expertise in biology and genetics.

    - This tutorial was written to be run from an AWS remote instance. You need an AWS account. Please see our [tutorial](../../Cloud-Platforms/Introduction_to_Amazon_Web_Services/introtoaws1.md) on setting up an AWS instance for help.

    - Basic shell scripting knowledge. Users must be comfortable with finding and opening a terminal window.

=== "Tutorial Resources"

    Vidlet: [BLAST tutorial walk-through](./BLAST2.md)
