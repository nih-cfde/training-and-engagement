---
layout: page
title: Terra Overview
hide:
  - toc
---

**An Introduction to Terra**

- what is Terra - cloud platform for data analysis
    - scale up
    - hosts existing datasets you might want to use
    - secure (even for human data) for sharing data/workflows/workspaces with collaborators
    - long-term persistence of data/workflow availability to other researchers

- what are the components of a workflow in Terra
    - a computer (GCP instance). For more about a GCP instance, see our [tutorial](../Introduction-to-GCP/index.md) on setting up a GCP instance.
    - software environment (docker)
    - a workflow of steps to do (WDL)
    - data (upload data to Terra workspace)



!!! bug "to figure out"

    how do we want ppl to follow this tutorial?

    do they create dockerfiles/WDLs that they actually upload to clouds? or just demonstrate how that happens?


Est. Time | Lesson name | Description
--- | --- | ---
30 mins  | [Setting up a Terra account](./1terra.md) |
30 mins  | [Navigating Terra workspaces](./2terra.md) |
30 mins  | [Running existing workflows](./3terra.md) |
1 hour  | [Custom workflow on Terra](./4terra.md) |

!!! note "Learning Objectives"

    The objectives of this tutorial are to:

    - set up Terra account

    - learn about the Terra platform interface

    - learn how to use existing workflows on Terra

    - learn how to upload a workflow to use on Terra

=== "Est. Cost"

    Testing and running workflows will incur some cost, but Google Cloud Platform gives new users a $300 free trial for 3 months.

    *Not sure how much it costs to test on GCP and/or how much actual run costs on Terra but it seems to be low*

=== "Prerequisites"

    For using Terra and its *existing* workflows:

    - GCP billing account (+ valid credit card). Please see our lesson for [setting up GCP accounts](../Introduction-to-GCP/index.md)

    For building *custom* workflows on Terra, you'll also need:

    - local terminal environment to test code - Mac and Linux terminal. Alternatively, all of this could be tested on a GCP instance (just incurs cost to run VMs, but if you are starting small it will not cost too much)
    - Dockerhub account
    - familiarity with shell command-line, docker, and Workflow Description Language computer programming




=== "Tutorial Resources"

    - [Terra support center](https://support.terra.bio/hc/en-us)
