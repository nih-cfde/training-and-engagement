---
layout: page
title: Terra Overview
hide:
  - toc
---

An Introduction to Terra
=========================

Terra is a cloud-based platform developed by the Broad Institute for bioinformatic analysis that aims to allow researchers with any level of computing skills to conduct large-scale bioinformatic analysis on the Google Cloud Platform.

Terra gives users both graphical and command line interface options, and also organizes documentation, data/metadata, workflows, analysis, and outputs all in one place. These shareable workspaces are important components for reproducible and collaborative research.

Est. Time | Lesson name | Description
--- | --- | ---
5 mins | [Introduction](0terra.md) | What is Terra and when would I use it?
30 mins  | [Setting up a Terra account](./1terra.md) | Connecting GCP and Terra accounts
30 mins  | [Navigating Terra workspaces](./2terra.md) | Intro to the Terra interface
30 mins  | [Cloud costs](./3terra.md) | How much does it cost to run an analysis?
1 hr | [Running existing workflows](./4terra.md) | How to run a workflow?
30 mins | [Custom workflow on Terra](./5terra.md) | Demo for building your own workflows

!!! note "Learning Objectives"

    The objectives of this tutorial are to:

    - learn about the Terra platform interface

    - set up a Terra account

    - learn how to use existing workflows on Terra

    - learn how to upload a workflow to use on Terra

=== "Est. Cost"

    Testing and running workflows will cost money, but the Google Cloud Platform gives new users a $300 free trial for 3 months.

    *reassess after deciding on demo analysis*

=== "Prerequisites"

    For using Terra and its **existing** workflows:

    - Google Cloud Platform (GCP) billing account (+ valid credit card). Please see our lesson for [setting up GCP accounts](../Introduction-to-GCP/index.md)

    For building **custom** workflows on Terra, you'll *also* need:

    - A terminal environment (i.e., Mac terminal or GCP virtual machine) to test workflow code. The GCP virtual machine will cost money however, so we recommend starting with small test workflows.
    - Dockerhub account (Dockstore is another option, but we will demonstrate with Dockerhub in this lesson)
    - Familiarity with shell command-line, docker, and Workflow Description Language computer programming are all needed to build custom workflows. The learning curve is steep, but there are also many resources available to get you started! We'll share a few in the lesson.

=== "Tutorial Resources"

    - [Terra support center](https://support.terra.bio/hc/en-us)
