---
layout: page
title: Rendering a GitHub website locally with Mkdocs
---

# Rendering a GitHub website locally with Mkdocs

This is a tutorial for building a website with mkdocs-material theme and hosting the website on readthedocs.com

*Python >2.7 is required. Consider installing inside a conda environment with python 3 or newer*

In the folder of your GitHub repo, you need a .yml file that specifies the layout of tabs on each webpage. All the actual webpages must be stored in a directory called `docs`. Finally, in the main folder you also need a text file that tells readthedocs to import the `material` theme. This file contains the following string: `mkdocs-material`.

## Install mkdocs using pip
If you are wondering what pip is and whether or not you must install it, please visit [this website](https://pip.pypa.io/en/stable/installing/).

To upgrade pip, run:

```
pip install --upgrade pip
```

### Install MkDocs

To install mkdocs, run this code:

```
pip install mkdocs
```
You will see some output. To check if mkdocs was successfully installed, run:

```
mkdocs --version
```
Your output should look something like this:

```
$ mkdocs --version
mkdocs, version 1.1.2
```

**Note**

If you are using Windows, some of the above commands may not work out-of-the-box.

A quick solution may be to preface every Python command with python -m like this:

```
python -m pip install mkdocs
python -m mkdocs
```

## Installing the material theme
You will use this theme to build your CFDE website. To import it, run:

```
pip install mkdocs-material
```

## Rendering locally

To render locally, run:

```
mkdocs serve
```

This takes ~30 seconds to run. Then copy the address shown in your command window ( it should look like "Serving on http://127.0.0.1:8000") and paste onto web browser.


## Hosting the website on readthedocs
*Note: you need admin and/or owner privileges to do this!*

Visit www.readthedocs.com

After log in, click `import project`. This should take you to a list of GitHub repos.

Click on the repo of interest.

Next go to `Project Settings` on the right nav-bar and click into `Advanced Settings`.


**Change the following setting:**

1. Set `Default version` to `latest`

1. Set `Default branch` to `name-of-your-theme`
    * You only need to do this step if you are changing the theme of the website. Readthedocs uses the master branch as the default branch to render the website. Since your new theme is on a different branch, it will not yet appear on your master branch. By changing the default branch in readthedocs, you can render your branch of choice.
    * If you are making minor changes to the theme (e.g. adding a file or fixing a typo), you may push changes to the [preview branch and use that link to render on a PR](https://cfde-training.readthedocs.io/en/latest/General-Tutorials/ProtectedBranch_HowTo/#preview-website-on-github-branch).

2. Scroll down to the `Documentation type` drop down menu and select `Mkdocs (Markdown)`

3. Click `Save` at the bottom of the page.

4. Now select the appropriate branch that you wish to build. In the case of major theme change, that branch would be `Latest`. And click the green button that says `Build`. This process takes a few mins.

5. Finally, make your website public like so:
    * Click on the `Versions` tab.
    * Click the `EDIT` button next to the `master` branch.
    * Find the setting called 'Privacy Level' and select 'Public' from the dropdown menu.

Enjoy your new website.
