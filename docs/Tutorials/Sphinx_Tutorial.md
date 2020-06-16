Sphinx for static websites editing
==================================

Install Sphinx for macOS
------------------------

Based on [this
tutorial](https://www.sphinx-doc.org/en/master/usage/installation.html)

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
\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~
\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~\~
