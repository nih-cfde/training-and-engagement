---
layout: page
title: Connect to your AWS instance on JupyterHub
---

# Connect to your AWS instance on JupyterHub

We've [installed some software for you on AWS](https://hackmd.io/rrjnYcZ3QemfuDpt82oomw?view), and you should receive a Web link to your AWS instance Web site in zoom chat (during today's demo, at least ;).

## Log in to JupyterHub

You should see:

![](https://i.imgur.com/h0mXlZy.png)

Use username `cfde`. (We'll give you the sekret password in Zoom chat.)

You should now be at a JupyerHub console.

![](https://i.imgur.com/W3o8BlN.png)

## Open Jupyter notebook

Go to the `demo/` folder and open `demo-walkthrough.ipynb` by clicking on it.

![](https://i.imgur.com/mWTz7DF.png)

You can [also see a static version of the walkthrough](https://github.com/ctb/2022-may-cfde-demo/blob/main/demo-walkthrough.ipynb).

## For exercise at end of Jupyter Notebook:

This will run through sample SRR5935743_1!

```
# run a sample against a database with sourmash gather
# take exactly 5.000 minutes
!sourmash gather SRR5935743_1.fastq.sig \
    gtdb-rs207.genomic-reps.dna.k31.zip \
      -o SRR5935743_1.gather.csv

# run the metacoder taxonomy plot
!./sourmash_gather_to_metacoder_plot.R \
    SRR5935743_1.gather.csv \
    SRR5935743_1.tax.png
```
