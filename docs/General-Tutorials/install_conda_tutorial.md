# Set up computing environment with conda on MacOS

Conda makes software installation and compute environment management easier. Conda can be installed via Miniconda or Anaconda. This tutorial is a walk-through with Miniconda.

> "Miniconda is a free minimal installer for conda. It is a small, bootstrap version of Anaconda that includes only conda, Python, the packages they depend on, and a small number of other useful packages, including pip, zlib and a few others. Use the conda install command to install 720+ additional conda packages from the Anaconda repository."

!!! tip

    Please refer to the [conda command cheatsheet](./basic_snakemake_tutorial/conda_cheatsheet.md) for commonly used conda commands!


## 1. Install conda

Follow the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to install Miniconda. Please use the 64-bit version when working with Python 3.x.

The steps for installing Miniconda for MacOS are [here](https://conda.io/projects/conda/en/latest/user-guide/install/macos.html).

### a. Download the [installer](https://docs.conda.io/en/latest/miniconda.html)
- Miniconda3 MacOSX 64-bit bash

### b. Verify your installer hashes
- go to the directory where you saved the installer file (e.g., Downloads/)
- open up Terminal and navigate to that directory
```
$ cd Downloads
$ shasum -a 256 Miniconda3-latest-MacOSX-x86_64.sh
```
- the output should look like this:
```
ccc1bded923a790cd61cd17c83c3dcc374dc0415cfa7fb1f71e6a2438236543d  Miniconda3-latest-MacOSX-x86_64.sh
```

### c. Install:
```
$ bash Miniconda3-latest-MacOSX-x86_64.sh
```
- you will be asked to review the license. Hit return key and scroll through.
- you will be asked if you accept the license terms. Type `yes`.
- Miniconda will be installed in the location printed on the screen. Hit return key to confirm.
- Miniconda will now be installed. This takes a few minutes to complete. The progress will be displayed in the Terminal window.
- you will be asked if you wish the installer to initialize Miniconda3 by running conda init. Type `yes`.
- for changes to take effect, close and re-open Terminal window.

### d. Verify that you can run conda
- now, when you re-open Terminal, the command prompt will start with `(base)`, indicating that you are in the base conda environment. It will look something like this:
```
(base) $
```
- try checking the version of your new conda installation:
```
(base) $ conda --version
```
- if you got something like `conda 4.8.3`, then you are ready for the next step!

### e. Configure conda
- conda uses channels to look for available software installations. These are some good channels to set up:
```
(base) $ conda config --add channels bioconda
(base) $ conda config --add channels conda-forge
```

## 2. Set up conda environment

There is always a `(base)` conda environment. You can then create new environments with different software set ups with the basic command:

```
(base) $ conda create -n <name of env>
```
This takes a few minutes (you'll see the message 'Solving environment'). Conda will then ask you to confirm the location of the new environment. Type `y`.

More options to customize the environment are documented under the help page for this command: `conda create -h`.

If you want to create an environment from a text file called 'environment.yml' that specifies the environment's requirements, the command would look like this:
`conda env create -n <new conda env name> -f environment.yml`. The `-f` flag specifies the .yml file that contains software requirements.

## 3. Activate conda environment

```
(base) $ conda activate <conda env name>`
```

Now, your command prompt starts with `(<conda env name>)`. If you named your new environment 'potato' it would look like this:
```
(potato) $
```

## 4. Take a look around your new conda environment!
This command shows you information about the conda environment you activated:
```
(potato) $ conda info
```

## 5. Install packages
The basic command for installing packages is:
```
(potato) $ conda install -y <software name>
```
It will ask if you want to install dependencies. Type `y`. This command will show a list of the software installed in this environment:
```
(potato) $ conda list -n <conda env name>
```

## 6. Leave conda environment

```
(potato) $ conda deactivate
```

You'll now be back in the `(base)` environment.
