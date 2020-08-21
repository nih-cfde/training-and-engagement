---
layout: page
title: BLAST Overview
---


# Running Command-Line BLAST

Learning Objectives

- Gain hands-on exposure to the linux command line
- Understand how data is turned into results by programs run at the command line

Tutorial Outline

- [BLAST Tutorial Video](../BLAST-Command-Line/BLAST1.md)
- [Install BLAST](../BLAST-Command-Line/BLAST2.md)
- [How to Run BLAST+](../BLAST-Command-Line/BLAST3.md)





## What is BLAST?

BLAST is the **B**asic **L**ocal **A**lignment **S**earch large sequence databases; It starts by finding small matches between the two sequences and extending those matches.  For in-depth information on how BLAST works and the different BLAST functionality, check out the resources page [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi).

BLAST can be helpful for identifying the source of a sequence or finding a similar sequence in another organism.  In this lesson, we will use BLAST to find zebrafish proteins that are similar to a small set of mouse proteins.

## Why use the command line?
BLAST has a very nice graphical interface for searching sequences in NCBI's database.
However, running BLAST through the commmand line has many benefits:


  * It's much easier to run many BLAST queries using the command line than the GUI
  * Running BLAST with the command line is reproducible and can be documented in a script
  * The results can be saved in a machine-readable format that can be analyzed later on
  * You can create your own databases to search rather than using NCBI's pre-built databases
  * It allows the queries to be automated
  * It allows you to use a remote computer to run the BLAST queries
  
Later on in the workshop we will talk more about these advantages and have a more in-depth explanation of the shell.

