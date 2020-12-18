---
layout: page
title: The Setup
---

Our goal is to analyze existing Kids First data in new ways,
or in new combinations, in order to improve medical outcomes. However,
before we can start using the data, we need to do a lot of set up.

Why so much setup?
------------------

There are many rules and regulations about who can
use human-derived medical data, and for what purposes. Even
"open" human data has controls and simply indicates that there are fewer barriers to seeing it. Data curators like
the Kids First DRC are obligated to enforce access controls and data use
rules to the best of their ability, and that requires end users like us
to have accounts that are tied to verified identities (e.g. ORCID or Google, eRA Commons).

Kids First uses their portal as a sort of catalog of their datasets, but
the data is stored in Cavatica, a platform for doing data analysis in the
cloud. We will need accounts for the Kids First
DRC Portal and Cavatica in order to do an analysis.


<!--Kids
First DRC maintains Whole Genome Sequences (WGS) and/or RNAseq data for
over 12,000 individuals. One type of file that stores genomic data like
this is called a bam file: a 	**B**inary sequence
**A**lignment **M**ap format file. We'll talk more
about file types later, but what is important here is that a bam file is
the smallest way to store alignment data.

A bam file for RNAseq from one
individual typically ranges from 15 to 30 *gigabytes*, while a WGS bam
file for one individual can be as many as 350GB. This does not include
all of the files that go with each bam file in order to make them
useable for analysis. As of early 2020, the Kids First overall dataset
is 1.31 *petabytes*. Since there is so much data, it needs to live in a
huge, dedicated compute space, and running an analysis generally
requires much more memory and storage than is available on an office
computer.
-->
