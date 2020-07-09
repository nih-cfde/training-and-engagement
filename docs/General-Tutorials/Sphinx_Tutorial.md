---
layout: page
title: Rendering a GitHub website locally with Sphinx
---

Sphinx for static websites editing
==================================

Install Sphinx for macOS
------------------------

Based on [this tutorial](https://www.sphinx-doc.org/en/master/usage/installation.html)

Create a new conda environment called `sphinx` that runs python 3.5 like
so:

    conda create --name sphinx python=3.5

Then activate the environment and installed sphinx with this code:

    conda activate sphinx
    conda install sphinx

Then navigate to your local GitHub folder that had the `make` file, and
run this:

*Note: this is the website I\'m working with: \`\`git clone
https://github.com/nih-cfde/training-and-engagement\`\`*

    make html

Then install dependencies based on error messages:

    pip install sphinx_bootstrap_theme
    pip install sphinx-markdown-tables
    pip install recommonmark

If you get error messages, you will need to run this code again after
all the installs:

    make html

If it worked, your html files can be found in a folder called:
`./build/html`

Click on individual .html pages to open them up in your web browser.

Learn how the various Sphinx files talk to each other in this [YouTube
video](https://www.youtube.com/watch?v=7adnbsj9A4w) by Paul Everitt


Install Sphinx for Windows
----------------------------

Based on this [tutorial](https://www.sphinx-doc.org/en/master/usage/installation.html)

**Step 1:** Install Python (*if necessary*)
Open `command prompt` for Windows

- To check if python 3 is installed on your computer, enter the command
```
python --version
```
If Python 3 is installed, you should see something similar to the following:
```
Python 3.7.4
```

- If Python is not installed, go to https://www.python.org/downloads/

Installing Python 3 enables you to install Sphinx with pip

**Step 2:** Install Sphinx - https://www.sphinx-doc.org/en/master/usage/quickstart.html

Open `command prompt` for Windows
```
 pip install -U sphinx
```
To check for sphinx, type 
```
sphinx-build --version
```
**Step 3:** To quickly build a local site 
```
sphinx-quickstart 
```
(It will generate a source directory with `conf.py` and a master document `index.rst`)

`index.rst` serves as a welcome page and to contain the root of the table of contents tree

**Step 4:** Open index.rst with your favorite text editor (Sublime 3, Atom, Notepad are some)

This allows you to get the basic template for the sphinx site. 

To create the html site, go to the root directory for your sphinx site. For me it is called
```
C:\Users\jsanc\sphinx-vidlet-hosting\docs
```
then type 
```
make html
``` 

Voila! You have a base template for a sphinx site.

