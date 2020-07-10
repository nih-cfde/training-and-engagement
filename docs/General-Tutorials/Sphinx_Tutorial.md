---
layout: page
title: Rendering a GitHub website locally with Sphinx
---

Sphinx for static websites editing
==================================

Install Sphinx for macOS
------------------------

Based on [this tutorial](https://www.sphinx-doc.org/en/master/usage/installation.html)

We will utilize conda package management system to install Sphinx. If you do not have conda installed, you can follow installation steps [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html).

Create a new conda environment called `sphinx` that runs python 3.5 like
so:

    conda create --name sphinx python=3.5

Then activate the environment and install Sphinx with this code:

    conda activate sphinx
    conda install sphinx

Build Sphinx site using template from GitHub
--------------------------------------------    

For this tutorial, we will use a template website hosted on GitHub, that was generated using Sphinx to modify and add content. First, let us create a local copy of the repo:

    git clone https://github.com/nih-cfde/Sphinx-demo.git

 Navigate to the newly created directory with the name of the repo.  To view the existing html template, navigate to `html` folder which is nested under `_build` folder, a subfolder of sphinx in the directory.
 Click on `index.html` to view the website on your browser.

 *Note that since GitHub does not allow for empty folders, some of Sphinx generated folders are missing from the repo. These can be directly added to the downloaded repo locally. The folders should be created within sphinx folder with exactly these name: `_static`, `_templates`, `docs`*.

 The configuration file `conf.py` is a python file that has details about the website rendering. The structure of the website is stored under `index.rst` which is in reStructured Text format.

 First we will modify existing text on the website. Open the `index.rst` file using any text editor. Remove the hyphen in the welcome header *Welcome to Sphinx-demo's documentation* and save.

 To add additional text, type under the welcome header after a single line spacing. *The line spacing is important part of syntax*. Here we added this statement and saved the `index.rst` file:

    This website documents building a static site using Sphinx
    as part of the Common Fund Data Ecosystem's (CFDE) training efforts.

 To render these changes, run the following code in the `sphinx` directory, which has the `Makefile`:

    make html

The newly rendered website with the added changes can be viewed by clicking `index.html` located in the `sphinx/_build/html` folder.

Adding files to Sphinx site
----------------------------  

Next, we will add some files to build the website. Let us create an empty file within the `docs` folder:

    touch installation.rst

Open the empty file using any editor. *Add '=' along the length of the heading on a new line to designate heading*. For example:

    Installing Sphinx
    =================

You can enter any text after a line spacing. Here, we added this:

    MacOS

    Sphinx can be installed using using `pip` or package management system like `conda`

Save and exit the file.
The `index.rst` has to be updated to render this file that we created. To do so, add `docs/installation.rst` after a single line spacing under `Contents:`.

Update the website by running the following code:

    make html

We can add multiple files to the website in similar manner.

Embed video on Sphinx site
---------------------------

In addition to text files, images and videos can also be embedded. Next, we will embed a talk by Paul Everitt exploring the feature set of Sphinx, hosted on [YouTube](https://www.youtube.com/watch?v=7adnbsj9A4w).

To embed YouTube video, click on `Share` and then `Embed` to obtain the `iframe` code. Accompanying text and the video are added above the `toctree` line in the `index.rst` file as shown:

    This talk by Paul Everitt highlights the various functionalities that can be achieved using Sphinx

    .. raw:: html

        <iframe width="560" height="315" src="https://www.youtube.com/embed/7adnbsj9A4w" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

Save, exit and compile using:

    make html

The rendered website will now have the embedded video on the landing page. Full set of Sphinx capabilities can be found at the official [Sphinx documentation site](https://www.sphinx-doc.org/en/master/).

Enjoy building websites!

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


