---
layout: page
title: MIME type Overview
hide:
  - toc
---

An Introduction to MIME types for File Formats
=================================================

A MIME type or media type is a form of identification for file formats and contents transmitted over the internet. It is useful to specify the data identification label of a file to allow software to properly interpret and render the data. This is especially important for Common Fund (CF) programs who may undertake data transfers over the internet and thus, have to ensure the data integrity along with data formats for a successful transfer. In this tutorial, we will describe how to determine MIME type for single and multiple files, and create custom MIME types specific to the file format.

Est. Time | Lesson name | Description
--- | --- | ---
5 mins | [Introduction](./Intro_MIME_type.md) | What is a MIME type ?
15 mins | [Example Data](./Example_data_files.md) | Obtain general or CF specific data files
5 mins | [file](./file.md) | MIME type using file
10 mins | [mimetype](./mimetype.md) | MIME type using mimetype
10 mins | [xdg-utils](./xdg-utils.md) | MIME type using xdg-utils
15 mins | [Siegfried](./siegfried.md) | MIME type using Siegfried
5 mins  | [Multi File MIME type ](./Multiple_file_MIME.md) | MIME type for multiple files
10 mins | [Unexpected Behavior](./Unexpected_behavior.md) | Understand and validate MIME type

!!! note "Learning Objectives"

    In this tutorial you will learn:

    - how to determine MIME type for single and multiple files

    - about different utilities for MIME type identification

    - how to create custom MIME types

    - how to validate MIME types for unknown file formats

=== "Prerequisites"

    This tutorial is written for a Unix/Linux compute environment (e.g. HPC, binder, cloud platforms like AWS, XSEDE etc). **The code included in the tutorial will not work on MacOS/X.** Some of the commands require super user privileges or permissions to use `sudo`.


=== "Tutorial Resources"

    Screencasts:

     - [file](./file.md)

     - [mimetype](./mimetype.md)

     - [xdg-utils](./xdg-utils.md)

     - [siegfried](./siegfried.md)
