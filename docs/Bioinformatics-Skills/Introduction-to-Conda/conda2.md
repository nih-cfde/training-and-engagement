# Let's get started!


- :computer: Open a web browser - e.g., Firefox, Chrome, or Safari.
- Open the binder by clicking this button: [![Binder](https://binder.pangeo.io/badge_logo.svg)](https://binder.pangeo.io/v2/gh/nih-cfde/training-rstudio-binder/conda-workshop-march2021?urlpath=rstudio)


Conda is already installed in the binder. We'll talk more about [setting conda up](./install_conda_tutorial.md) on your local system later in the lesson!




Set up Rstudio panels!

### Initialize conda

:keyboard: Copy/paste commands into the terminal OR run the commands from the `workshop_commands.sh` file in the binder.

Installer sets up two things: Conda and the root environment. The root environment contains the selected python version and some basic packages.

![](https://i.imgur.com/NM4wtre.png)
Image credit: [Gergely Szerovay](https://www.freecodecamp.org/news/why-you-need-python-environments-and-how-to-manage-them-with-conda-85f155f4353c/)

Setup the conda installer and initialize the settings:

```bash=1
conda init
```

We will shorten command prompt to `$`:
```bash=2
echo "PS1='\w $ '" >> .bashrc
```

Re-start terminal for the changes to take effect (type `exit` and then open a new terminal).

:::success
:heavy_check_mark: Put up a :hand: on Zoom when you've got a new terminal window and see this at the prompt:
`~ $`
:::

---

### Conda channels: Searching for software

The channels are places that conda looks for packages. The default channels after conda installation is set to Anaconda Inc's channels (Conda's Developer).

```bash=3
conda config --show channels
conda list # get list of packages in base environment
```

Channels exist in a hierarchial order. By default:
:::info
Channel priority > package version > package build number
:::


![](https://i.imgur.com/GE8nY0a.png)
Image credit: [Gergely Szerovay](https://www.freecodecamp.org/news/why-you-need-python-environments-and-how-to-manage-them-with-conda-85f155f4353c/)

:::info
 See [Resources](https://hackmd.io/1ivIgqVYSCupC21MuAW9FA?view#Resources) for more info on channels!
:::

We will update the channel list order and add `bioconda` since we are using bioinformatic tools today. **The order of the channels matters!**

First, add the `defaults` channel:

```bash=5
conda config --add channels defaults
conda config --get channels
```

Then add the `bioconda` channel:
```bash=7
conda config --add channels bioconda
conda config --get channels
```

Lastly, add the `conda-forge` channel to move it to top of the list, following [Bioconda's recommended channel order](https://bioconda.github.io/user/install.html#set-up-channels):

```bash=9
conda config --add channels conda-forge
conda config --get channels
```

With this configuration, conda will search for packages in this order: 1) `conda-forge`, 2) `bioconda`, and 3) `defaults`

:::info
Another way to add channels is:

```
conda config --prepend channels bioconda
```
:::

:question: Questions?

---

**Instructor switch!**

### Install Software

We will install FastQC which is a software tool that provides a simple way to run quality control checks on raw sequencing data.


Search for software (fastqc):

```bash=1
conda search fastqc
```

Create conda environment and install fastqc:

```bash=2
conda create -y --name fqc fastqc
```

Activate environment:

```bash=3
conda activate fqc
```
Check fastqc version:

```bash=4
fastqc --version
```
:::success
:heavy_check_mark: Put up a :hand: on Zoom if your command prompt shows the name of the environment in parentheses:
`(fqc) ~ $`
:::

:::info
To go back to `(base) ~ $` environment:
```
conda deactivate
```
:::

---

High-throughput sequencing data quality control steps often involve fastqc and trimmomatic. Trimmomatic is useful for read trimming (i.e., adapters). There are multiple ways we could create a conda environment that contains both software programs:

#### Method 1: install software in existing environment

We could add `trimmomatic` to the `fqc` environment:
```bash=5
conda install -y trimmomatic=0.36
conda list # check installed software
```

We can specify the exact software version :point_up_2:
The default is to install the most current version, but sometimes your workflow may depend on a different version.

---

#### Method 2: install both software during environment creation

When you switch conda environments, conda changes the PATH (and other environment variables) so it searches for software packages in different places.

Let's check the PATH for method 1:
```
echo $PATH
```
You should see that the first element in the PATH changes each time you switch environments!

```bash=
conda deactivate
conda create -y --name fqc_trim fastqc trimmomatic=0.36
conda activate fqc_trim
conda list # check installed software
echo $PATH # path for method 2
```
---

The following methods use an external file to specify the packages to install.

#### Method 3: specify software to install with a YAML file
Often, it's easier to create environments and install software using a YAML file (YAML is a file format) that specifies all the software to be installed. For our example, we are using a file called `test.yml`. Let's start back in the `(base)` environment.

```bash=1
conda deactivate
```

The `test.yml` file contains the following in YAML format:

```yaml
name: qc_yaml #this specifies environment name
channels:
    - conda-forge
    - bioconda
    - defaults
dependencies:
    - fastqc
    - trimmomatic=0.36
```

Create the environment - note the difference in conda syntax:
```bash=2
conda env create -f test.yml #since environment name specified in yml file, we do not need to use -n flag here
conda activate qc_yaml
conda list  # check installed software
```

:::info
See [Resources](https://hackmd.io/RweX2WZ7RGGtfLfCDuEubQ?view#How-do-I-recreate-an-environment-with-exactly-the-same-versions-of-all-software-packages) for a 4th method of exporting exact environments.

Software can also be installed by specifying the channel with `-c` flag i.e.:
```
conda install -c conda-forge -c bioconda sourmash
```
:::

---

### Managing Environments

At this point, we have several conda environments! To see a list, there are 2 commands (they do the same thing!):
```bash
conda env list
```

or

```bash
conda info --envs
```

!!! warning
Generally, you want to avoid installing too many software packages in one environment. It takes longer for conda to resolve compatible software versions for an environment the more software you install.

For this reason, in practice, people often manage software for their workflows with multiple conda environments.



---

### Running FastQC in the fqc Environment

If not already done, activate one of the environments we created, e.g.,:

```bash
conda activate fqc
```

Let's make sure the software was installed correctly:
```bash
fastqc --help
```

:::success
Output should look like:
```
FastQC - A high throughput sequence QC analysis tool

SYNOPSIS

        fastqc seqfile1 seqfile2 .. seqfileN

    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
           [-c contaminant file] seqfile1 .. seqfileN
...
```
:::

Download data (a yeast sequence file):

```bash=
curl -L https://osf.io/5daup/download -o ERR458493.fastq.gz
```

Check out the data:
```bash=2
gunzip -c ERR458493.fastq.gz | wc -l
```


What does the fastq file look like?
```bash=3
gunzip -c ERR458493.fastq.gz | head
```

Run FastQC!
```bash=4
fastqc ERR458493.fastq.gz
```


## Exercise 1 :writing_hand:

- **First**, install the `blast` software with version 2.9.0 using conda.

:::spoiler Hint!
:bulb: You can use one of these approaches to install `blast`:
   - install the software in an existing env using `conda install -y <name of the software>`
   - create a new env using `conda create -y --name <name of env> <software to install>`
:::

- **Second**, try this example BLAST analysis in the conda environment. In this example, we are comparing the sequence similarity of a mouse protein sequence to a reference zebra fish sequence - as you might imagine, they are not that similar!
But for today, this exercise will demonstrate running a quick analysis in a conda environment and bonus points if you find out how similar/dissimilar they are! (More details on BLAST and what each step is for [here](https://training.nih-cfde.org/en/latest/Bioinformatics-Skills/Command-Line-BLAST/BLAST4/)).
Run each line of code below in the terminal:

    - Make a directory for the exercise files:
    ```bash=
    mkdir exercise1
    cd exercise1
    ```
    - Download with `curl` command and unzip data files:
    ```bash=3
    curl -o mouse.1.protein.faa.gz -L https://osf.io/v6j9x/download
    curl -o zebrafish.1.protein.faa.gz -L https://osf.io/68mgf/download
    gunzip *.faa.gz
    ```
    - Subset the data for a test run:
    ```bash=6
    head -n 11 mouse.1.protein.faa > mm-first.faa
    ```
    - Format zebra fish sequence as the blast database to search against:
    ```bash=7
    makeblastdb -in zebrafish.1.protein.faa -dbtype prot
    ```
    - Run a protein blast search with `blastp`!
    ```bash=8
    blastp -query mm-first.faa -db zebrafish.1.protein.faa -out mm-first.x.zebrafish.txt -outfmt 6
    ```

What does the output look like?

:::info
Note that if you `conda deactivate`, you can still access the input/intermediate/output files from the BLAST analysis. They are not 'stuck' inside the conda environment!
:::


---

## Exercise 2 :writing_hand:

Conda allows you to revert to a previous version of your software using the `--revision` flag:

Usage:
```
# list all revisions
conda list --revisions

# revert to previous state
conda install --revision <number>

# for example:
conda install --revision 1
```

Earlier, we installed an older version of trimmomatic (0.36). Try updating it to the most recent version and then revert back to the old version.

:::spoiler Hint!
:bulb: You can do this exercise in any of the conda environments we created earlier with trimmomatic. You can update software with `conda update <software name>`
:::

---

**Tidying up**

To remove old conda environments:
```
conda env remove --name <conda env name>
```

To remove software:
```
conda remove <software name>
```

:::warning
Be sure to save any work/notes you took in the binder to your computer. Any new files/changes are not available when the binder session is closed!

For example, select a file, click "More", click "Export":
![](https://i.imgur.com/js8eX0l.png)
:::
