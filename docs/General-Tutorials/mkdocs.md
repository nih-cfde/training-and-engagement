---
layout: page
title: Rendering a GitHub website locally with MkDocs
---

# Rendering a GitHub website locally with MkDocs

MkDocs is a static site generator built for project documentation. It comes with easy to use and customizable themes and features. The documentation files are written in Markdown and configured using a single configuration file in YAML format. The generated site can be hosted on any hosting platform including [GitHub pages](https://pages.github.com), [Amazon S3](https://docs.aws.amazon.com/AmazonS3/latest/dev/WebsiteHosting.html) or [readthedocs](https://readthedocs.com/).  

!!! note "Learning Objectives"

    - Install MkDocs
    - Add website content
    - Change default theme
    - Host the website on [readthedocs.com](https://readthedocs.com/)  

=== "Est. Time"

    35 mins

=== "Prerequisites"

    - Access to MacOS, Windows or Linux system
    - To install MkDocs the system should either have Python > v3.5 or suitable package manager like `conda` installed.
    - Have [`git`](https://git-scm.com/) installed on your computer
    - Basic command line skills
    - Admin and/or owner access to [readthedocs.com](https://readthedocs.com/) is required for hosting the website.  

=== "Tutorial Resources"    

    [Demo GitHub repository with MkDocs Template](https://github.com/nih-cfde/MkDocs-demo)   

## Step 1: Install MkDocs

MkDocs can be installed via a package manager or manually using `pip` the Python package manager. In this tutorial, installation using `pip` and `conda` are detailed. Choose either one installation route to begin.

=== "pip install"

    `pip` is python package manager. To use `pip` as an installer you will need [Python](https://www.python.org/downloads/) installed on your system. For more details on `pip` installation, please visit [the official website](https://pip.pypa.io/en/stable/installing/). The instructions to set up MkDocs using `pip` follow the [tutorial from MkDocs](https://www.mkdocs.org/#installation).

    !!! note "MkDocs prerequisite"
        Python v3.5 or higher is required. Consider installing inside a `conda` environment with python v3.5 or newer.

    To upgrade pip, run:

    ```
    pip install --upgrade pip
    ```

    Then install MkDocs:

    ```
    pip install mkdocs
    ```

    To check if MkDocs was successfully installed, run `mkdocs` with `--version` flag:

    ```
    $ mkdocs --version
    mkdocs, version 1.1.2
    ```

    !!! note "Modification for Windows OS"
        Some of the installation commands may not compile correctly on Windows OS. Running the python module as a script might fix it. To do so, add `-m` flag to python commands:  

         ```
         python -m pip install mkdocs
         python -m mkdocs
         ```

=== "conda install"

    We can also utilize conda package management system to install MkDocs. If you do not have conda installed, you can follow installation steps for [MacOS](./install_conda_tutorial.md),[Windows](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html) or [Linux](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) operating systems.

    Create a new conda environment called `mkdocs` that runs the latest version of python 3:

    ```
    conda create --name mkdocs python=3
    ```
    Then activate the environment:

    ```
    conda activate mkdocs
    ```

    MkDocs is hosted on `conda-forge` channel. We first add this channel to our conda configuration:

    ```
    conda config --add channels conda-forge
    ```

    To ensure we are installing the latest version of MkDocs, we can search all available versions:

    ```
    conda search mkdocs --channel conda-forge
    ```

    Install the latest version which was `v1.1.2` at the time of this tutorial:

    ```
    conda install mkdocs=1.1.2
    ```

## Step 2: Build MkDocs site using template from GitHub

For this tutorial, we will use a template website hosted on GitHub, that was generated using MkDocs to modify and add content. First, let us create a local copy of the repo:

```
git clone https://github.com/nih-cfde/MkDocs-demo.git
```

Navigate to the newly created directory with the name of the repo. The `.yml` file in `yaml` format specifies the layout of tabs on each webpage. All the actual webpages are stored in the `docs` directory and currently it contains an "index.md" file.

We can render the website locally.

=== "Input"
    ```
    mkdocs serve
    ```

=== "Expected Output"
    ```
    # Generated output on a MacOS

    (mkdocs) scanchi@ mkdocs-demo $ mkdocs serve
    INFO    -  Building documentation...
    INFO    -  Cleaning site directory
    INFO    -  Documentation built in 0.06 seconds
    [I 200713 13:58:47 server:334] Serving on http://127.0.0.1:8000
    INFO    -  Serving on http://127.0.0.1:8000
    [I 200713 13:58:47 handlers:62] Start watching changes
    INFO    -  Start watching changes
    [I 200713 13:58:47 handlers:64] Start detecting changes
    INFO    -  Start detecting changes
    ```

Copy and paste the server address to a web browser to render the site. When you are done checking the local version, ++ctrl+c++ to close the server.

## Step 3: Add content on MkDocs site

One key feature of the dev-server that MkDocs offers is the auto-reloading when any change is detected. Open the "mkdocs.yml" file using any text editor and change site name to `MkDocs Trial` and save. You will notice the website will automatically reflect the change in title.

Next we will contents for `About` section of the website. Create a blank markdown file in `docs` folder and name it "about.md". Add some sample contents to this file and save:

```
# About Mkdocs

This site was generated using MkDocs v1.1.2
```

We have to add navigation information to the configuration file which will dictate the order, title and nesting of the additional pages. The updated "mkdocs.yml" file will contain:

```
site_name: MkDocs Trial
nav:
    - Home: index.md
    - About: about.md
```

These changes are reflected in the website with the `Home`, `About` icons on the left along `Search`, `Previous` and `Next` on the right side on the top navigation bar. The navigation bar additions and the search features are integrated in MkDocs without requiring additional configuration on the user end.

## Step 4: Change default theme

MkDocs comes installed with two themes: `Bootstrap` and `readthedocs`. To change the theme to `readthedocs`, add this line in the "mkdocs.yml" file:

```
theme: readthedocs
```

To recreate the look of [CFDE training website](https://cfde-training-and-engagement.readthedocs-hosted.com/en/latest/), we will install an external MkDocs theme called [Material for MkDocs](https://github.com/squidfunk/mkdocs-material). List of all available external MkDocs themes can be found at [MkDocs wiki page](https://github.com/mkdocs/mkdocs/wiki/MkDocs-Themes). Choose either `pip` or `conda` to install  the material theme.

=== "pip install"

    ```
    pip install mkdocs-material
    ```

=== "conda install"

    ```
    conda install mkdocs-material
    ```

We can either choose to install a specific version for this theme or by default install the latest available version. The website should automatically update with the new theme.

Options to deploy the generated website are elaborated on the [MkDocs official website](https://www.mkdocs.org/user-guide/deploying-your-docs/).

For the website to be hosted on [readthedocs.com](https://readthedocs.com/) a text file with instructions to import the `material` theme needs to the added to the main folder of the website. This file contains the following string: `mkdocs-material`.

## Step 5: Host the website on Read the Docs

!!! note "Read the Docs requirement"
    You need admin and/or owner privileges to host the website.

Visit [readthedocs.com](https://readthedocs.com/)

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

!!! note "Key Points"

    You now have a workflow to customize your static website, generate content and host on [readthedocs.com](https://readthedocs.com/). Enjoy your new website!
