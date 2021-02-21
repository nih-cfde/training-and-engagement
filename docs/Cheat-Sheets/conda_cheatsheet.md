# Conda Command Cheat Sheet

Commonly used conda commands:

conda | Description
--- | ---
`conda create -n <conda env name>` or `conda create -n <conda env name> <software name>` | Create a new conda environment. You can include other flags to customize the environment more and install software in the conda environment. For example, `conda create -n fastqc_env fastqc` will install the FastQC program into a conda environment called `fastqc_env`.
`conda env create -n <conda env name> -f <yaml file>` | Create a new conda environment using specifications from a yaml file. For example, `conda env create -n test -f environment.yml`.
`conda install -y <software name>` | Install software in conda environment
`conda activate <conda env name>` | Activate conda environment
`conda deactivate` | Deactivate conda environment
`conda info --envs` or `conda env list` | Both commands list conda environments, `*` will be next to the environment you are currently in
`conda list -n <conda env name>` | List software installed in this conda environment. Or simply, `conda list`.
`conda info` | Information about your conda environment
`conda search <software name>` | Search for available software versions
`conda env remove --name <conda env name>` | Remove a conda environment


Download the official [conda cheat sheet](https://docs.conda.io/projects/conda/en/latest/user-guide/cheatsheet.html) for more commands
