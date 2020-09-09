## Tutorial components

### Tutorial audience

Who is our target audience?
- clinicians, biomedical data scientists 
- NIH Common Fund staff
- NIH Common Fund Data Ecosystem Coordinating Center staff

In general, tutorials should avoid assuming the user's experience level, unless specified as a prerequisite.

For a thorough guide to tutorial style best practices, please read [DigitalOcean's style guide section](https://www.digitalocean.com/community/tutorials/digitalocean-s-technical-writing-guidelines#style).

### Tutorial content

Tutorials should consist primarily of original content. If lesson material is adapted from other sources, please attribute work appropriately (e.g., "This material was adapted from ANGUS (link to source material)."). If you want to include original video or animation content from other sources, check that:

1) the material has a license that allows reuse, 
2) if there is a license, the material is properly referenced in the tutorial,
3) where possible, we tell the original authors we are using their material and appreciate their work! and,
4) the addition of these materials are supplements to the tutorial and not the primary content.

### Landing page

Multi-page tutorials should have a summary page of the tutorial prerequisites and contents. Landing pages are not necessary for 1-page tutorials, however the same information should be provided in the tutorial's introduction section instead. 

Tutorial landing page components:
- Title
- Short introduction
- Learning Objectives
- Prerequisites
- Table of tutorial pages
- Tutorial resources

See the [landing page template](./tutorial_template_docs/TutorialLandingPg_template.md) for a page outline.

*Title*
Titles should be short and include the goal of the tutorial, e.g., `How to launch an AWS instance`, `How to create a website with Mkdocs`

*Short introduction*
A brief description of what the tutorial is about (expanded in the Introduction section, see [Tutorial structure](#tutorial-structure))

*Learning Objectives

The objectives should address:
- what are the key points users should learn?
- what are the expected outcomes of the tutorial (e.g., new skills learned)?

They are formatted with a `note` admonition box:

```
!!! note "Learning Objectives"
```

*Prerequisite*

Clearly state the operating system(s) that will work for the tutorial, any required installation or set up steps that are not documented in the tutorial (e.g., you can link to existing set up tutorials instead). While we recommend writing tutorials without making assumptions about experience level, if particular computational, bioinformatics, and/or biological knowledge are needed, state them clearly as prerequisites.

*Table*

A Markdown formatted table for each page of the tutorial. The table columns should be named as follows:

Est. Time | Lesson name | Description
--- | --- | ---

Est. Time: Provide estimated completion times for each page (e.g., 30 secs, 10 mins, 1 hour)
Lesson name: The title of each page, with a hyperlink to the page
Description: Formatted as a short phrase or question about the primary learning goal of each section

*Tutorial resources*
Add in hyperlinks to any tutorial reference material (e.g., cheatsheets, scripts, example data, vidlets, etc.)


### Tutorial structure

The following sections should be included in each tutorial:

- Title
- Introduction
- Set up
- Tutorial steps 
- Conclusion

See the [tutorial main page template](./tutorial_template_docs/TutorialMainPg_template.md) for a page outline.

*Title*

For tutorials that have a landing page, the title will refer to the primary content of each page. For tutorials that do not have a separate landing page, titles should be short and include the goal of the tutorial, e.g., `How to launch an AWS instance`, `How to create a website with Mkdocs`.

*Introduction*

The introduction should briefly address:
- what is the tutorial about (general topic, e.g., cloud computing with AWS, workflows with Snakemake)?
- what will the user do in this specific tutorial (e.g., create an account on Kids First, create a Manhattan plot)?

*Set up*

Provide instructions on computer set up, software installations, and/or file downloads (may include linking to existing tutorials).

*Steps*

Tutorial steps should be written as `Step n: Action`, e.g., `Step 1: Downloading reference file`. 

Include links to accompanying files for download, like scripts or input files for data downloads

For command-line software installations, provide a test command to check that the installation was successful. This could be checking for the software version or help page:
```
# input
conda --version
# output
conda 4.8.3

# input
conda --help
# truncated output
usage: conda [-h] [-V] command ...

conda is a tool for managing and deploying applications, environments and packages.

Options:
```

*Conclusion*

Briefly sum up the tutorial key points as a bullet point list (it could be a restatement of the learning objectives if any were stated in the Introduction). Use the admonition format:

```
!!! note "Key Points"
    
    - point1
    - point2
```


**Saving files**
- save markdown files in `./docs/`
- save images in `./docs/images/` as .jpeg or .png files
- there are several ways to save files that accompany tutorials:
    - save reference material in `./docs/Resources/` (e.g., code syntax cheatsheets)
    - for example data files that users should access with `wget` or `curl`, upload files to the [CFDE Training Data Files osf project](https://osf.io/c8txv/#show_login). Login to the osf project is restricted to CFDE members. Create a new folder and upload your file(s).
    - for small files that users should download from the website (e.g., template script files), format the file name as a hyperlink in the tutorial. If the file has a file extension, a download option is available when users click the link.
    - for scripts/example scripts associated with a specific tutorial, save files in the tutorial folder.
    - save binder files in a new Github repo (binders are made for entire repos, and we do not want a binder of the entire training-and-engagement repo!)

For future discussion:
- pending need: save scripts in a general `./docs/tutorial_scripts/` folder 

### Additional/optional user-friendly components
- walk-through vidlets
- walk-through screencasts
- binder set-up so users can follow along without worrying about computer set up or software installations
- set up practice Github repository

 
