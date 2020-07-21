---
layout: page
title: Rendering a GitHub website locally with Sphinx
---

Sphinx for static websites editing
==================================

Sphinx is python based tool that enables generation of documentation. Originally intended for python projects, it has versatility to be adapted for any number of projects that generate plain text.

Install Sphinx for macOS
------------------------

Based on [this tutorial](https://www.sphinx-doc.org/en/master/usage/installation.html).

We will utilize conda package management system to install Sphinx. If you do not have conda installed, you can follow installation steps [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html).

Create a new conda environment called `sphinx` that runs python 3.5 like
so:

```
conda create --name sphinx python=3.5
```

Then activate the environment and install Sphinx with this code:

```
conda activate sphinx
conda install sphinx
```

Install Sphinx for Windows
----------------------------

Based on this [tutorial](https://www.sphinx-doc.org/en/master/usage/installation.html).

Install Python (*if necessary*). *Note you can also use conda to install both python and Sphinx in Windows*.

Open `command prompt` for Windows.
To check if python 3 is installed on your computer, enter the command:

```
$ python --version
Python 3.7.4
```

If Python is not installed, you can download from [here](https://www.python.org/downloads/).

To install [Sphinx](https://www.sphinx-doc.org/en/master/usage/quickstart.html), we will use `pip` which is a python package installer.

Open `command prompt` for Windows and type the following code:

```
pip install -U sphinx
```

The argument `U` is to enable upgrade of the packages.

To check if Sphinx was installed, type on `command prompt`:

```
sphinx-build --version
```

Build Sphinx site using template from GitHub
--------------------------------------------    

For this tutorial, we will use a template website hosted on GitHub, that was generated using Sphinx to modify and add content. First, we will create a local copy of the repo:

```
git clone https://github.com/nih-cfde/Sphinx-demo.git
```

Navigate to the newly created directory with the name of the repo.  

!!! note "Add folders"
    Since GitHub does not allow for empty folders, some of Sphinx generated folders are missing from the repo. These can be directly added to the downloaded repo locally. The folders should be created within sphinx folder with exactly these names: `_static`, `_templates`, `docs`.

This is the overall structure of the sphinx directory after the necessary folders have been created:

```
├── Makefile
├── _build
│   ├── doctrees
│   │   ├── environment.pickle
│   │   └── index.doctree
│   └── html
│       ├── _sources
│       ├── _static
│       ├── genindex.html
│       ├── index.html
│       ├── objects.inv
│       ├── search.html
│       └── searchindex.js
├── _static
├── _templates
├── conf.py
├── docs
├── index.rst
└── make.bat
```

To view the existing html template, navigate to `html` folder nested under `_build` folder, a subfolder of sphinx in the directory. Click on `index.html` to view the website on your browser.

The configuration file `conf.py` is a python file that has details about the website rendering. The structure of the website is stored under `index.rst` , which is in reStructured Text format and serves as the welcome page for the website.

First we will modify existing text on the website. Open the `index.rst` file using any text editor. Remove the hyphen in the welcome header *Welcome to Sphinx-demo's documentation* and save.

To add additional text, type under the welcome header after a single line spacing. *The line spacing is important part of syntax*. The following sample text was added  to `index.rst` and saved:

```
This website documents building a static site using Sphinx
as part of the Common Fund Data Ecosystem's (CFDE) training efforts.
```

To render these changes, run the following code in the `sphinx` directory, which has the `Makefile`:

```
make html
```

The newly rendered website with the added changes can be viewed by clicking `index.html` located in the `sphinx/_build/html` folder.

Adding files to Sphinx site
----------------------------  

Next, we will add some files to build the website. Let us create an empty file within the `docs` folder:

```
touch installation.rst
```

Open the empty file using any editor. *Add '=' along the length of the heading on a new line to designate heading*. For example:

```
Installing Sphinx
=================
```

You can enter any text after a line spacing. Here, we added this:

```
MacOS

Sphinx can be installed using using `pip` or package management system like conda
```

Save and exit the file.
The `index.rst` has to be updated to render this file that we created. To do so, add `docs/installation.rst` after a single line spacing under `Contents:`.

Update the website by running the following code:

```
make html
```

We can add multiple files to the website in similar manner.

Embed video on Sphinx site
---------------------------

In addition to text files, images and videos can also be embedded. Next, we will embed a talk by Paul Everitt exploring the feature set of Sphinx, hosted on [YouTube](https://www.youtube.com/watch?v=7adnbsj9A4w).

To embed YouTube video, click on `Share` and then `Embed` to obtain the `iframe` code. Accompanying text and the video are added above the `toctree` line in the `index.rst` file as shown:

```
This talk by Paul Everitt highlights the various functionalities that can be achieved using Sphinx

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/7adnbsj9A4w" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```

Save, exit and compile using:

```
make html
```

The rendered website will now have the embedded video on the landing page. Full set of Sphinx capabilities can be found at the official [Sphinx documentation site](https://www.sphinx-doc.org/en/master/).

Build local Sphinx site
------------------------

A previously created website template was used in this tutorial. However, you can also build a local version of the Sphinx site.

To do so, on the command line, type:

```
sphinx-quickstart
```

This will bring up a set of queries regarding the setting up of the source directory. You can choose the default values.

This will generate `conf.py` and `index.rst` along with necessary folders which will be empty.

To generate the html site, compile in the same directory that contains the `Makefile`

```
make html
```

The `index.html` will be accessible in the `html` folder nested under `docs` folder. This will generate the basic template using Sphinx. You can add contents, change configuration of the website and customize following the steps discussed in this tutorial.

Enjoy building websites!
