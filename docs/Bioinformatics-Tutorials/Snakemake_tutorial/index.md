---
layout: page
title: Snakemake Overview
---

**An introduction to Snakemake for workflow management**

Workflow management systems help to automate analyses and make them easier to maintain, reproduce, and share with others. In this tutorial, we will walk through the basic steps for creating a [variant calling](https://www.ebi.ac.uk/training-beta/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/) workflow with the Snakemake workflow management system.

This may not be the variant calling workflow you would necessarily use in practice, but it serves as a good example for teaching Snakemake. Many people do indeed use `samtools`, but for particularly big or complex genomes, guidelines provided by  [GATK](https://gatk.broadinstitute.org/hc/en-us) would serve best. Additionally, various parameters associated with mapping, visualization etc may require tuning.

!!! goal

    The objectives of this tutorial are to:
    
    - learn how to write basic workflows with Snakemake rules
    
    - learn variable substitution for Snakemake rules
    
    - learn wildcard matching for Snakemake rules
    
    - understand why workflow systems can help you do your computing more easily
    
These materials were was adapted from DIB lab course materials [here](https://github.com/ngs-docs/2020-GGG298) and [here](https://github.com/ngs-docs/2020-GGG201b-lab).

!!! note "Prerequisites"
    
    This tutorial is written for a Unix or Linux compute environment (e.g., MacOS, Linux-based HPC, pre-configured binder). It assumes basic knowledge of navigating, editing files, and executing scripts from the command line. Some knowledge of Python is useful, but not required.    

Est. Time | Lesson name | Description
--- | --- | ---
5 mins | [Introduction](./snakemake_0.md) | What is a workflow? <br />What is Snakemake?
30 mins | [Setup](./snakemake_1.md) | Set up tutorial computing environment
30 mins | [The Snakefile](./snakemake_2.md) | What is the Snakefile? <br />What are Snakemake rules?
45 mins | [Decorating the Snakefile](./snakemake_3.md) | How do the rules link together?
30 mins | [Continue Decorating](./snakemake_4.md) | More on rule linking
20 mins | [Run Through the Workflow](./snakemake_5.md) | Final Snakemake rule <br />What are the results?

Tutorial resources:

- Videos: 

    - Part 1: [Introducing the Snakefile and Snakemake](./snakemake_2.md)

    - Part 2: [Continue decorating the Snakefile](./snakemake_4.md)

    - Part 3: [More decorating and running through the entire Snakemake workflow](./snakemake_4.md)

- [Example Snakefile](./example_snakefile.md)

- Cheatsheets:

    - [bash and nano](./bash_cheatsheet.md)
    
    - [conda](./conda_cheatsheet.md)
  
    - [Snakemake](./snakemake_cheatsheet.md)
