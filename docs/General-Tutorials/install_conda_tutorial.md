# Set up computing environment with conda on MacOS

Conda makes software installation and compute environment management easier by making sure that all the software you are using for a particular project works together. Conda can be installed via Miniconda (a smaller more efficient package) or Anaconda (the full installation of conda). This tutorial is a walk-through with Miniconda. From the conda website:

> "Miniconda is a free minimal installer for conda. It is a small, bootstrap version of Anaconda that includes only conda, Python, the packages they depend on, and a small number of other useful packages, including pip, zlib and a few others. Use the conda install command to install 720+ additional conda packages from the Anaconda repository."

Est. Time | Lesson name 
--- | --- 
20 mins | Set up conda environment 

!!! note "Learning Objectives"

    Learn how to install conda and set up a conda environment.

=== "Prerequisites"

    This tutorial is written specifically for installing the MacOS version of Miniconda. 
    
=== "Tutorial Resources"

    Please refer to the [conda command cheatsheet](../Resources/conda_cheatsheet.md) for commonly used conda commands!


Follow the Miniconda [installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). Please use the 64-bit version when working with Python 3.x.

There are specific steps to install Miniconda for [MacOS](https://conda.io/projects/conda/en/latest/user-guide/install/macos.html).

### Step 1: Download the [installer](https://docs.conda.io/en/latest/miniconda.html)

Select "Miniconda3 MacOSX 64-bit bash"

### Step 2: Verify your installer hashes

Go to the directory where you saved the installer file (e.g., "Downloads/"). Open up Terminal and navigate to that directory:

```
cd Downloads
```

=== "Input"

    ```
    shasum -a 256 Miniconda3-latest-MacOSX-x86_64.sh
    ```

=== "Expected Output"
    
    ```
    ccc1bded923a790cd61cd17c83c3dcc374dc0415cfa7fb1f71e6a2438236543d  Miniconda3-latest-MacOSX-x86_64.sh
    ```

### Step 3: Install Miniconda

```
bash Miniconda3-latest-MacOSX-x86_64.sh
```

- you will be asked to review the license. Hit `RETURN` key and scroll through.

- you will be asked if you accept the license terms. Type `yes`.

- Miniconda will be installed in the location printed on the screen. Hit `RETURN` key to confirm.

- Miniconda will now be installed. This takes a few minutes to complete. The progress will be displayed in the Terminal window.

- you will be asked if you wish the installer to initialize Miniconda3 by running conda init. Type `yes`.

- for changes to take effect, `exit` and re-open Terminal window.

### Step 4: Verify that you can run conda

Now, when you re-open Terminal, the command prompt will start with `(base)`, indicating that you are in the base conda environment. It will look something like this:

```
(base) $
```

Check the version of your new conda installation:

=== "Input"

    ```
    conda --version
    ```
    
=== "Expected Output"

    If you got a conda version, then you are ready for the next step!
    
    ```
    conda 4.8.3
    ```

### Step 5: Configure conda

Conda uses channels to look for available software installations. These are some good channels to set up:

```
conda config --add channels bioconda
```

```
conda config --add channels conda-forge
```

### Step 6: Set up conda environment

There is always a `(base)` conda environment. You can then create new environments with different software set ups with the basic command:

```
conda create -n <name of env>
```

This takes a few minutes (you'll see the message "Solving environment"). Conda will then ask you to confirm the location of the new environment. Type `y`.

More options to customize the environment are documented under the help page for this command: `conda create -h`.

If you want to create an environment from a text file called "environment.yml" that specifies the environment's requirements, the command would look like this:
`conda env create -n <new conda env name> -f environment.yml`. The `-f` flag specifies the `.yml` file that contains software requirements.

### Step 7: Activate conda environment

```
conda activate <conda env name>
```

Now, your command prompt starts with `(<conda env name>)`. If you named your new environment "potato" it would look like this:

```
(potato) $
```

### Step 8: Take a look around your new conda environment!
This command shows you information about the conda environment you activated:

```
conda info
```

### Step 9: Install packages
The basic command for installing packages is:

```
conda install -y <software name>
```

It will ask if you want to install dependencies. Type `y`. This command will show a list of the software installed in this environment:

```
conda list -n <conda env name>
```

### Step 10: Leave conda environment

```
conda deactivate
```

You'll now be back in the `(base)` environment.

!!! note "Key Points"

    Now you should have a working conda installation that you can use to create custom conda environments!

