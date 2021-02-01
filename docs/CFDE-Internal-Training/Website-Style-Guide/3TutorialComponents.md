# Tutorial components

## Tutorial audience

Who is our target audience?

- clinicians, biomedical data scientists

- NIH Common Fund staff

- NIH Common Fund Data Ecosystem Coordinating Center staff

In general, tutorials should avoid assuming the user's experience level, unless specified as a prerequisite.

For a thorough guide to tutorial style best practices, please read [DigitalOcean's style guide section](https://www.digitalocean.com/community/tutorials/digitalocean-s-technical-writing-guidelines#style).

## Tutorial content

Tutorials should consist primarily of original content. If lesson material is adapted from other sources, please attribute work appropriately (e.g., "This material was adapted from ANGUS (link to source material)."). If you want to include original video or animation content from other sources, check that:

1) the material has a license that allows reuse,

2) if there is a license, the material is properly referenced in the tutorial,

3) where possible, we tell the original authors we are using their material and appreciate their work! and,

4) the addition of these materials are supplements to the tutorial and not the primary content.

## Tutorial structure

All tutorials should begin with landing page information. For longer tutorials that are split over multiple pages, start the tutorial steps on a new page (more details below). For one-page tutorials, the landing "page" information may all be on the same page as the tutorial steps.

!!! tip

    See the [markdown tutorial template](https://github.com/nih-cfde/training-and-engagement/blob/dev/docs/CFDE-Internal-Training/Website-Style-Guide/tutorial_template_docs/TutorialTemplate.md) for a page outline. You can copy/paste the markdown to start new tutorials.

### Tutorial landing page components

In this order:

- Title

- Brief description

- Table of contents

- "Learning Objectives"

- "Prerequisites"

- "Tutorial Resources" (optional)

*Title*

Titles should be short and include the goal of the tutorial, e.g., `How to launch an AWS instance`, `How to create a website with MkDocs`

*Brief description*

A brief description of what the tutorial is about (may be expanded in the Introduction section detailed below)

*Table of contents (for MULTI-page tutorials)*

A Markdown formatted table for each page of the tutorial. The table columns should be named as follows:

Est. Time | Lesson name | Description
--- | --- | ---
Est. Time: Provide estimated completion times for each page (e.g., 30 secs, 10 mins, 1 hour) | Lesson name: The title of each page, with a hyperlink to the page | Description: Formatted as a short phrase or question about the primary learning goal of each section

*"Learning Objectives"*

The objectives should address:

- what are the key points users should learn?

- what are the expected outcomes of the tutorial (e.g., new skills learned)?

They are formatted with a `note` admonition box:

```
!!! note "Learning Objectives"
```

*Est. Time (for 1-page tutorials)*

For **1-page** tutorials, include the "Est. Time" as a tabbed box rather than a table and put it AFTER the "Learning Objectives" so that the time, prereqs, and tutorial resources are grouped together (see Conda installation tutorial for an example):

```
=== "Est. Time"

    10 mins
```

Where applicable (e.g., cloud computing tutorials), include an "Est. Cost" section as a tabbed box:

```
=== "Est. Cost"

    < $1.00 to run through entire tutorial.

```

*"Prerequisites"*

Clearly state the operating system(s) that will work for the tutorial, any required installation or set up steps that are not documented in the tutorial (e.g., you can link to existing set up tutorials instead). While we recommend writing tutorials without making assumptions about experience level, if particular computational, bioinformatics, and/or biological knowledge are needed, state them clearly as prerequisites.

They are formatted with a tabbed box:

```
=== "Prerequisites"
```

*"Tutorial Resources" (optional)*

If applicable, add in hyperlinks to any tutorial reference material (e.g., cheat sheets, scripts, example data, vidlets, etc.).

They are formatted with a tabbed box:

```
=== "Tutorial Resources"
```

### Sections after landing page

**The following sections (in this order) should be included in each tutorial, following the landing page, either on the same page or on a new page for longer tutorials:**

- Title (optional)

- "Introduction" (optional)

- "Set Up" (optional)

- Tutorial steps

- Conclusion (optional)

*Title (optional)*

For tutorials that start new pages after their landing page, the title will refer to the primary content of each page, e.g., `How to launch an AWS instance`, `How to create a website with MkDocs`.

*"Introduction" (optional)*

If more introductory material is needed beyond the brief description on the landing page, please add an introduction section that addresses:

- what is the tutorial about (general topic, e.g., cloud computing with AWS, workflows with Snakemake)?

- what will the user do in this specific tutorial (e.g., create an account on Kids First, create a Manhattan plot)?

*"Set Up" (optional)*

If more set up material is needed beyond the prerequisites on the landing page, please provide additional instructions in a Set Up section (e.g., on computer set up, software installations, and/or specific tutorial file downloads (may include linking to existing tutorials)).

*Steps*

Tutorial steps should be written as `Step n: Action`, e.g., `Step 1: Downloading reference file`.

Include links to accompanying files for download, like scripts or input files for data downloads

For command-line software installations, **provide a test command** to check that the installation was successful. This could be checking for the software version or help page:

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

*Conclusion (optional)*

If applicable, briefly sum up the tutorial key points as a bullet point list (it could be a recap of the learning objectives). In some cases, it may make sense to have a "Conclusion" section and in others, a "Key Points" admonition box is sufficient.

Use the admonition `note` format:

```
!!! note "Key Points"

    - point1
    - point2
```

### Saving files

- save markdown files in `./docs/` in the appropriate tutorial folder ("Bioinformatics-Skills", "CFDE-Internal-Training", "Web-Development")

- save images as .jpeg or .png files. Images should be saved in a subfolder within the relevant tutorial folder. For example, the Snakemake tutorial is located in `./docs/Bioinformatics-Skills/Snakemake` and the images used in the tutorial are located in `./docs/Bioinformatics-Skills/Snakemake/images-snakemake`.

- there are several ways to save files that accompany tutorials:

    - save reference material in `./docs/Cheat-Sheets/` (e.g., general code syntax cheat sheets)

    - for example data files that users should access with `wget` or `curl`, upload files to the [CFDE Training Data Files osf project](https://osf.io/c8txv/#show_login). Login to the osf project is restricted to CFDE members. Create a new folder and upload your file(s).

    - for small files that users should download from the website (e.g., template script files), format the file name as a hyperlink in the tutorial. Give the file a file extension (e.g., "testfiles.txt"), so a download option is available when users click the link.

    - for scripts/example scripts associated with a specific tutorial, save files in the tutorial folder.

    - for vidlets, see the [recording style guide](./4RecordingStyleGuide.md)

    - save binder files in a new Github repo (binders are made from entire repos, and we do not want a binder of the entire training-and-engagement repo!).

## Additional/optional user-friendly components

Check out the next pages for guidelines on the following:

- walk-through [vidlets and screencasts](./4RecordingStyleGuide.md)

- [binder set-up](./5PangeoBinderGuide.md) so users can follow along without worrying about computer set up or software installations

- set up [practice Github repository](./6PracticeGithubRepos.md)
