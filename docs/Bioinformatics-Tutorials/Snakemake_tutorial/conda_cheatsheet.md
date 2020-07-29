# conda command cheatsheet

conda | Description
--- | ---
`conda create -n <conda env name>` | create a new conda environment. You can include other flags to customize the environment more.
`conda activate <conda env name>` | activate conda environment
`conda install -y <software name>` | install software in conda environment
`conda deactivate` | deactivate conda environment
`conda info --envs` | list conda environments, `*` will be next to the environment you are currently in
`conda list -n <conda env name>` | list software installed in this conda environment. Or simply, `conda list`.
`conda info` | information about your conda environment
`conda env remove --name <conda env name>` | remove a conda environment
