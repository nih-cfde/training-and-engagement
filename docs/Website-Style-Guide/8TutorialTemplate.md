---
layout: page
title: Template
---

YAML header:

    Overview/landing pages must start with the yaml file header, e.g.,:

    ```
    ---
    layout: page
    title: <tutorial name> overview
    hide:
       - toc
    ---
    ```
    This sets the title of the previous and upcoming pages on the Mkdocs footer bar.

    To ensure the contents occupy the full width of the page use:

    ```
    hide:
      - toc
    ```

    in the yaml format on top of the page.

Landing page title
==================

- Use the `===` under the titles for landing page.

Add a brief description of what the tutorial is about. If tutorial material was built from other sources, mention that here.

For **multi-page** tutorials only, add estimated time per section, lesson name, and description in a Markdown table:

Est. Time | Lesson name | Description
--- | --- | ---
x mins | [Page title](path/to/page) | Formatted as a short phrase or question about the primary learning goal of each section
y mins | [Page title](path/to/page) | Shorter tutorials may have only one entry

!!! note "Learning Objectives"

    The objective(s) of this tutorial are to:

    - objective 1

    - objective 2

For **1-page** tutorials only, add a tabbed box for the "Est. Time" (it should be grouped with the prereqs and tutorial resources tabs, and as necessary, the estimated cost tab):

=== "Est. Time"

    x mins

=== "Est. Cost"

    $1.00

=== "Prerequisites"

    - Operating system(s)
    - Required installations
    - Setup steps needed to complete the tutorial

=== "Tutorial Resources"

    - list out resources with bullet points and include a hyperlink to the resource (e.g., vidlets, screencasts, example files, cheat sheets)


**For longer tutorials that are split over multiple pages, start the next sections on a NEW page. Use the `#` header style on new pages.**

### Introduction (optional)

If more introductory material is needed beyond the brief description on the landing page, please add an introduction section.

### Set Up (optional)

If more set up material is needed beyond the prerequisites section, please provide additional instructions in a Set Up section (e.g., on computer set up, software installations, and/or specific tutorial file downloads (may include linking to existing tutorials)).

### Step 1: do this

### Step 2: do that

### Step 3: change this

### Conclusion (optional)

If applicable, briefly sum up the tutorial key points as a bullet point list (it could be a recap of the learning objectives). In some cases, it may make sense to have a "Conclusion" section and in others, a "Key Points" admonition box is sufficient.

Use the admonition format:

!!! note "Key Points"

    - point1
    - point2
