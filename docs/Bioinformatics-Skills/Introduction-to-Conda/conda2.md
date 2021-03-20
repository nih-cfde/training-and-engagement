# Let's get started!


- :computer: Open a web browser - e.g., Firefox, Chrome, or Safari.
- Open the binder by clicking this button: [![Binder](https://binder.pangeo.io/badge_logo.svg)](https://binder.pangeo.io/v2/gh/nih-cfde/training-rstudio-binder/conda-workshop-march2021?urlpath=rstudio)


!!! info

    For this lesson, we are using 3 of the Rstudio panels: Source panel to run conda commands, Terminal panel to execute code, and File panel to view input/output files.

Conda is already installed in the binder. We'll talk more about [setting conda up](./install_conda_tutorial.md) on your local system later in the lesson!


!!! warning

    What happens if I get a 502, 503, or 504 error from the binder?

    Try clicking on the launch button again. The binder or internet connection may have timed out.




### Initialize conda

:keyboard: Copy/paste commands into the terminal OR run the commands from the "workshop_commands.sh" file in the binder (in File Rstudio panel).

Installer sets up two things: Conda and the root environment. The root environment contains the selected python version and some basic packages.

![](./conda-imgs/conda-init.png 'conda installer')

Image credit: [Gergely Szerovay](https://www.freecodecamp.org/news/why-you-need-python-environments-and-how-to-manage-them-with-conda-85f155f4353c/)

Setup the conda installer and initialize the settings:

```
conda init
```

We will shorten command prompt to `$`:
```
echo "PS1='\w $ '" >> .bashrc
```

Re-start terminal for the changes to take effect (type `exit` and then open a new terminal).




There is always a `(base)` conda environment.



### Conda channels: Searching for software

The channels are places that conda looks for packages. The default channels after conda installation is set to Anaconda Inc's channels (Conda's Developer).

```
conda config --show channels
conda list # get list of packages in base environment
```

Channels exist in a hierarchial order. By default:

Channel priority > package version > package build number

![](./conda-imgs/conda-channel.png 'conda channels')

Image credit: [Gergely Szerovay](https://www.freecodecamp.org/news/why-you-need-python-environments-and-how-to-manage-them-with-conda-85f155f4353c/)

!!! info

    Commonly used channels:

    - In absence of other channels, conda [searches the `defaults` repository](https://docs.anaconda.com/anaconda/user-guide/tasks/using-repositories/) which consists of ten official repositories
    - `conda-forge` and `bioconda` are channels that contain community contributed software
    - `Bioconda` specializes in bioinformatics software (*supports only 64-bit Linux and Mac OS*)
        - [package list](https://anaconda.org/bioconda/repo)
    - `conda-forge` contains many dependency packages
        - [package list](https://anaconda.org/conda-forge/repo)
    - You can even install R packages with conda!

We will update the channel list order and add `bioconda` since we are using bioinformatic tools today. **The order of the channels matters!**

First, add the `defaults` channel:

```
conda config --add channels defaults
conda config --get channels
```

Then add the `bioconda` channel:
```
conda config --add channels bioconda
conda config --get channels
```

Lastly, add the `conda-forge` channel to move it to top of the list, following [Bioconda's recommended channel order](https://bioconda.github.io/user/install.html#set-up-channels). This is because many packages on bioconda rely on dependencies that are available on conda-forge, so we want conda to search for those dependencies before trying to install any bioinformatics software.

```
conda config --add channels conda-forge
conda config --get channels
```

With this configuration, conda will search for packages in this order: 1) `conda-forge`, 2) `bioconda`, and 3) `defaults`

!!! info

    Another way to add channels is:

    ```
    conda config --prepend channels bioconda
    ```

### Install Software

We will install [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) which is a software tool that provides a simple way to run quality control checks on raw sequencing data.


Search for software (fastqc):

```
conda search fastqc
```

Create conda environment and install FastQC. This takes a few minutes (you'll see the message "Solving environment").

```
conda create -y --name fqc fastqc
```


More options to customize the environment are documented under the help page for this command: `conda create -h`.

Activate environment:

```
conda activate fqc
```

This command shows you information about the activated conda environment:

```
conda info
```

Check fastqc version:

```
fastqc --version
```

!!! info

    To go back to `(base) ~ $` environment:
    ```
    conda deactivate
    ```

High-throughput sequencing data quality control steps often involve FastQC and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). Trimmomatic is useful for read trimming (i.e., adapters). There are multiple ways we could create a conda environment that contains both software programs:

#### Method 1: install software in existing environment

We could add `trimmomatic` to the `fqc` environment:
```
conda install -y trimmomatic=0.36
conda list # check installed software
```

We can specify the exact software version :point_up_2:
The default is to install the most current version, but sometimes your workflow may depend on a different version.

!!! info

    Software can also be installed by specifying the channel with `-c` flag i.e.:
    ```
    conda install -c conda-forge -c bioconda sourmash
    ```

#### Method 2: install both software during environment creation

When you switch conda environments, conda changes the PATH (and other environment variables) so it searches for software packages in different folders.

Let's check the PATH for method 1:
```
echo $PATH
```
You should see that the first element in the PATH changes each time you switch environments!

```
conda deactivate
conda create -y --name fqc_trim fastqc trimmomatic=0.36
conda activate fqc_trim
# check installed software
conda list
# path for method 2
echo $PATH
```

The following methods use an external file to specify the packages to install:

#### Method 3: specify software to install with a YAML file
Often, it's easier to create environments and install software using a YAML file that specifies all the software to be installed. For our example, we are using a file called `test.yml`. Let's start back in the `(base)` environment.

```
conda deactivate
```

The `test.yml` file contains the following in YAML format:

=== "YAML"
```
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
```
conda env create -f test.yml #since environment name specified in yml file, we do not need to use -n flag here
conda activate qc_yaml
conda list  # check installed software
```

!!! info

    [YAML](https://en.wikipedia.org/wiki/YAML) is a file format that is easy for both computers and humans to read. YAML file extensions are `.yml` and these files can be generated in any text editor.

#### Method 4: Install exact environment

For this approach, we export a list of the exact software package versions installed in a given environment and use it to set up new environments. This set up method won't install the latest version of a given program, for example, but it will replicate the exact environment set up you exported from.

```
conda activate fqc
conda list --export > packages.txt
conda deactivate
```

Two options -

1) install the exact package list into an existing environment:

```
conda install --file=packages.txt
```

OR

2) set up a new environment with the exact package list:
```
conda env create --name qc_file --file packages.txt
```

### Managing Environments

At this point, we have several conda environments! To see a list, there are 2 commands (they do the same thing!):
```
conda env list
```

OR

```
conda info --envs
```

Note that the current environment you're in is marked with an asterisk `*`.

!!! warning

    Generally, you want to avoid installing too many software packages in one environment. It takes longer for conda to resolve compatible software versions for an environment the more software you install.

    For this reason, in practice, people often manage software for their workflows with multiple conda environments.

### Running FastQC in a conda environment

Let's run a small analysis with FastQC in the `fqc` environment we created above.

If not already done, activate one of the environments we created, e.g.,:

```
conda activate fqc
```

Let's make sure the software was installed correctly:

=== "Input"

    ```
    fastqc --help
    ```

=== "Expected Output"

    Output should look like:
    ```
    FastQC - A high throughput sequence QC analysis tool

    SYNOPSIS

        fastqc seqfile1 seqfile2 .. seqfileN

        fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
           [-c contaminant file] seqfile1 .. seqfileN
    ...
    ```


Download data (a yeast sequence file):

```
curl -L https://osf.io/5daup/download -o ERR458493.fastq.gz
```

Check out the data:

=== "Input"

    The `gunzip -c` command allows us to see the unzipped version of the file without actually unzipping it (you can verify this by checking the file extension after running this command!). The `|` is called a pipe and it takes the output of the `gunzip -c` command and hands it to the `wc` word count command. The `-l` flag tells `wc` we want to count the number of lines in the file.

    ```
    gunzip -c ERR458493.fastq.gz | wc -l
    ```

=== "Expected Output"

    There should be 4,375,828 lines in the file.

What does the fastq file look like?

=== "Input"

    Here again, we use the `gunzip -c` and pipe the output to the `head` command to show the first 10 lines of the file:

    ```
    gunzip -c ERR458493.fastq.gz | head
    ```

=== "Expected Output"

    The beginning of the fastq format sequence file should look like this, where the 1st line is the sequence read ID (starts with `@`), the 2nd line is the DNA sequence, the 3rd is sequence separator `+`, and the 4th is the Phred quality score associated with each base pair in ASCII format.

    ```
    @ERR458493.1 DHKW5DQ1:219:D0PT7ACXX:1:1101:1724:2080/1
    CGCAAGACAAGGCCCAAACGAGAGATTGAGCCCAATCGGCAGTGTAGTGAA
    +
    B@@FFFFFHHHGHJJJJJJIJJGIGIIIGI9DGGIIIEIGIIFHHGGHJIB
    @ERR458493.2 DHKW5DQ1:219:D0PT7ACXX:1:1101:2179:2231/1
    ACTAATCATCAACAAAACAATGCAATTCAAGACCATCGTCGCTGCCTTCGC
    +
    B@=DDFFFHHHHHJJJJIJJJJJJIJJJJJJJJJJJJJJJJJJJJIJJJJI
    @ERR458493.3 DHKW5DQ1:219:D0PT7ACXX:1:1101:2428:2116/1
    CTCAAAACGCCTACTTGAAGGCTTCTGGTGCTTTCACCGGTGAAAACTCCG
    ```

Run FastQC!

=== "Input"

    ```
    fastqc ERR458493.fastq.gz
    ```

=== "Expected Output"

    On the terminal screen, FastQC prints analysis progress:
    ```
    Started analysis of ERR458493.fastq.gz
    Approx 5% complete for ERR458493.fastq.gz
    Approx 10% complete for ERR458493.fastq.gz
    Approx 15% complete for ERR458493.fastq.gz
    Approx 20% complete for ERR458493.fastq.gz
    Approx 25% complete for ERR458493.fastq.gz
    Approx 30% complete for ERR458493.fastq.gz
    Approx 35% complete for ERR458493.fastq.gz
    ...
    Analysis complete for ERR458493.fastq.gz
    ```

    The final output file is called "ERR458493_fastqc.html".

    You can click on the `.html` file in the File panel to open it in a web browser. This is the quality check report for our yeast sequence file.

## Let's practice!

=== "Exercise 1"

    - **First**, install the `blast` software with version 2.9.0 using conda.

    - **Second**, try this example BLAST analysis in the conda environment. In this example, we are comparing the sequence similarity of a mouse protein sequence to a reference zebra fish sequence - as you might imagine, they are not that similar!
    But for today, this exercise will demonstrate running a quick analysis in a conda environment and bonus points if you find out how similar/dissimilar they are! (More details on BLAST and what each step is for [here](../Command-Line-BLAST/BLAST4.md)). Run each line of code below in the terminal:

        - Make a directory for the exercise files:
        ```
        mkdir exercise1
        cd exercise1
        ```

        - Download with `curl` command and unzip data files:
        ```
        curl -o mouse.1.protein.faa.gz -L https://osf.io/v6j9x/download
        curl -o zebrafish.1.protein.faa.gz -L https://osf.io/68mgf/download
        gunzip *.faa.gz
        ```

        - Subset the data for a test run:
        ```
        head -n 11 mouse.1.protein.faa > mm-first.faa
        ```

        - Format zebra fish sequence as the blast database to search against:
        ```
        makeblastdb -in zebrafish.1.protein.faa -dbtype prot
        ```

        - Run a protein blast search with `blastp`!
        ```
        blastp -query mm-first.faa -db zebrafish.1.protein.faa -out myoutput.txt -outfmt 6
        ```

    What does the output (`myoutput.txt`) look like?

=== "Hint"

    :bulb: You can use one of these approaches to install `blast`:

    - install the software in an existing env using `conda install -y <name of the software>`

    - create a new env using `conda create -y --name <name of env> <software to install>`

=== "Answer"

    There are many ways to solve this - here is one approach!

    Input:
    ```
    # search for blast versions
    conda search blast
    # create an environment and install blast version 2.9.0
    conda create -y -n exercise1 blast=2.9.0
    # activate environment
    conda activate exercise1

    # make a folder for exercise 1 files
    mkdir exercise1
    # go to the exercise1 folder
    cd exercise1

    # download input sequence files
    curl -o mouse.1.protein.faa.gz -L https://osf.io/v6j9x/download
    curl -o zebrafish.1.protein.faa.gz -L https://osf.io/68mgf/download
    # unzip the gzip files
    gunzip *.faa.gz

    # make a small example subset of the query sequence file
    head -n 11 mouse.1.protein.faa > mm-first.faa
    # format the reference database
    makeblastdb -in zebrafish.1.protein.faa -dbtype prot
    # run a protein BLAST search
    blastp -query mm-first.faa -db zebrafish.1.protein.faa -out myoutput.txt -outfmt 6
    ```

    Output looks like this:
    ```
    YP_220550.1	NP_059331.1	69.010	313	97	0	4	316	10	322	1.24e-150	426
    YP_220551.1	NP_059332.1	44.509	346	188	3	1	344	1	344	8.62e-92	279
    YP_220551.1	NP_059341.1	24.540	163	112	3	112	263	231	393	5.15e-06	49.7
    YP_220551.1	NP_059340.1	26.804	97	65	2	98	188	200	296	0.10	35.8
    ```

    The 3rd column has the percent of matching amino acids between the query mouse protein sequence file and the zebra fish reference database. Not surprisingly, :mouse: and :fish: are not very similar! For this particular mouse sequence, it's only 69.01% similar to the zebra fish reference.

    Note that if you `conda deactivate`, you can still access the input/intermediate/output files from the BLAST analysis. They are not 'stuck' inside the conda environment!

---

=== "Exercise 2"

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

=== "Hint"

    :bulb: You can do this exercise in any of the conda environments we created earlier with trimmomatic. You can update software with `conda update <software name>`

=== "Answer"

    There are many ways to solve this - here is one approach!

    ```
    # activate fqc environment
    conda activate fqc
    # check version of trimmomatic
    conda list
    # update to latest version, should be 0.39 or higher
    conda update trimmomatic
    # look at revision list
    conda list --revisions
    # install the revision to revert to, here we want to revert to (rev 1) which had 0.36
    conda install --revision 1
    # now it's back to 0.36!
    conda list
    ```

    Note that conda told us that it was upgrading to 0.39 and downgrading to 0.36 - useful outputs!



### Tidying up

To remove old conda environments:
```
conda env remove --name <conda env name>
```

To remove software:
```
conda remove <software name>
```

!!! warning

    Be sure to save any work/notes you took in the binder to your computer. Any new files/changes are not available when the binder session is closed!

    For example, select a file, click "More", click "Export":
    ![](./conda-imgs/binder-save-files.png "save binder files")
