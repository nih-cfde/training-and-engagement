---
layout: page
title: Overview
hide:
  - toc
---

An Introduction to UNIX for Cloud Computing
============================================

This tutorial introduces the UNIX command line. It is designed for scientists and clinicians who need to use cloud-based and remote computers for basic and biomedical research.

!!! note Learning Objectives

	The objectives of this tutorial are to:

	* Understand basic UNIX command structure 
	* Navigate through hierarchical directory structures
	* Read, write, create, copy, move, and remove files
	* Understand wildcards and regular expression 
	* Redirect outputs and write for loops to build reproducible workflo	

| Lesson name | Description | Commands
| --- | --- | --- |
| [The shell and terminal](./unix_1.md) | What is UNIX? <br /> Why should I use it? | `PS1='$ '`, `clear` |
| [Navigating files and directories](./unix_2.md) | Where am I? <br /> What files are here? | `pwd`, `cd`, `ls` | 
| [Reading large files](./unix_3.md) | How do I read compressed files? <br />  Are my files of good quality?  | `head`, `tail`, `cat`, `less`, `gunzip`, `zip` |
| [Creating and modifying files](./unix_4.md) | How do I combine commands to build workflows? | `cp`, `mv`, `mkdir`, `rm`, `curl` |
| [Finding things](./unix_5.md) | How do I save my results? <br /> What commands did I type? | `grep`, `find` |
| [Redirection and for loops](./unix_6.md) | How do I combine commands to build workflows?   <br /> Do my results match published? | `wc`, `echo`, `for`
| [Exploring .csv files](./unix_7.md) | Is my gene of interest present this dataset?  <br /> How often? | `cut` | 
| [Concluding thoughts](./unix_8.md) | How do I save my history? <br /> What resources do you recommend? | `history` | 


=== "Tutorial Resources"

Please refer to our [UNIX Cheatsheet](https://training.nih-cfde.org/en/latest/General-Tools/Cheat-Sheets/bash_cheatsheet/) for a list of commonly used commands.

The website [Explain Shell](https://explainshell.com/) searches the help text and prints the information for a command line entered. 

The lesson materials were adapted from the UC Davis Data Lab's [Intro to Cloud Computing](https://ngs-docs.github.io/2021-august-remote-computing/) workshop, Data Carpentry's [Introduction to the Command Line for Genomics](https://datacarpentry.org/shell-genomics/) lesson, and the Lab for Data Intensive Biology's [Advanced Beginner/Intermediate Shell](https://dib-training.readthedocs.io/en/pub/2016-01-13-adv-beg-shell.html) workshop. 
