# What else can we use Conda for?

- **Reproducing analyses**
    - **Why**: removes the guesswork on what version of software to use or how to install it!

    :::spoiler
    For example, it's a lot easier to replicate a software environment from a publication or share software versions used in your own analysis using conda (particularly method 3 & 4 mentioned above ways to create conda environments from YAML or exact package list files).
    :::

- **Packaging software**
    - **Why**: to share your cool software with the world!

    :::spoiler
    The gist is that you write a recipe that has all the specs about your software that is submitted to a channel, e.g., conda-forge or bioconda. Once correctly formatted and tested (with continuous integration automation) so it works on different operating systems, it's added to the channel and the world can use and install your software with conda!
    - https://python-packaging-tutorial.readthedocs.io/en/latest/conda.html
    :::

- **Running parts of analysis workflows in their own environment**
    - **Why**: avoid software version conflicts, easier to keep track of software/versions used for a specific analysis, easier for others to reproduce

    :::spoiler
    This is a very helpful feature for workflows. For example, in [Snakemake](https://training.nih-cfde.org/en/latest/Bioinformatics-Skills/Snakemake/) workflows, you can specify that each step (called a "rule") is executed in an isolated conda environment by adding a `conda:` directive:

    ```python
    rule fastqc_raw:
        input: "rnaseq/raw_data/{sample}.fq.gz"
        output: "rnaseq/raw_data/fastqc/{sample}_fastqc.html"
        params:
             outdir="rnaseq/raw_data/fastqc"
        conda: "rnaseq-env.yml"
        shell:
              """
              fastqc {input} --outdir {params.outdir}
              """
    ```
    :::

- **Creating binders**
    - **Why**: teaching tool for software or analysis demos, share reproducible analysis

    :::spoiler
    Examples:

    - the Rstudio binder we're using today was created with https://binder.pangeo.io/
    - [metENP binder](https://github.com/metabolomicsworkbench/MetENP): tool for metabolite enrichment analysis and their associated enriched pathways
    :::

## Wrapping up

There are many actions you can perform with conda environments. Today we have mentioned these!
- init
- config
- search
- create
- activate/deactivate
- list
- remove
- update
- revert
- export




## Resources

- [CFDE training webite](https://training.nih-cfde.org/en/latest/)
- Common conda commands:
    - [CFDE conda cheat sheets](https://training.nih-cfde.org/en/latest/Cheat-Sheets/conda_cheatsheet/)
    - [official conda cheat sheet](https://docs.conda.io/projects/conda/en/latest/user-guide/cheatsheet.html)
- [Workshop additional resources/FAQ](https://hackmd.io/RweX2WZ7RGGtfLfCDuEubQ)
- [CFDE events & workshops](https://www.nih-cfde.org/events/)
- [DIB Lab video about conda](https://www.youtube.com/watch?v=Ef1QwhELuMs)
