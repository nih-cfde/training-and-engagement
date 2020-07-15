# Decorating the Snakefile

But what if we want to run all the commands at once? It gets tedious to run each command individually, and we can do that already without snakemake!

By defining the inputs and outputs for each rule's command/commands, snakemake can figure out how the rules are linked together. The rule structure will now look something like this:
```
rule <rule name>:
    input:
        # input file names must be enclosed in quotes
        # multiple inputs should be separated by commas
        # the new line for each input is optional
        <input file>,
        <input file 2>,
        <input file 3>
    output:
        # output file names must be enclosed in quotes
        # multiple outputs should be separated by commas
        <output file>,
        <output file 2>
    shell:
        # for multi-line commands
        # commands must be enclosed in triple quotes
        """
        <command 1>
        <command 2>
        """
```

Here, snakemake interprets the `input:` and `output:` sections as Python code, and the `shell:` section as the bash code that gets run on the command line.

### Add input and output files

Let's start with a clean slate. Delete any output files you created in the sections above, such that you only have the Snakefile in your directory: `rm <file name>`. Be careful with this command - it deletes files forever!

**Adding outputs:**

The output of the `download_data` rule is 'SRR2584857_1.fastq.gz'. Add this to the rule, note that the output file must be in quotes `""`:
```
rule download_data:
    output: "SRR2584857_1.fastq.gz"
    shell:
        "wget https://osf.io/4rdza/download -O SRR2584857_1.fastq.gz"
```

Try: Run the `download_data` rule twice.

Hey, nothing happened!
![](../../images/snakemake_nothingtobedone.jpeg)

Try removing the file: `rm SRR2584857_1.fastq.gz`. Now run the rule again.

It ran! This is because snakemake is smart enough to know that an output file already exists and doesn't need to be re-created.

**Adding inputs:**

To the `download_genome` rule, add: `output: "ecoli-rel606.fa.gz"`

and to the `uncompress_genome` rule, add an input and output:
```
rule uncompress_genome:
    input: "ecoli-rel606.fa.gz"
    output: "ecoli-rel606.fa"
    shell:
        "gunzip ecoli-rel606.fa.gz"
```
What does this do?

This tells snakemake that `uncompress_genome` depends on having the input file `ecoli-rel606.fa.gz` in this directory, and that `download_genome` produces it. What's extra cool is that snakemake will now automatically run the rule that outputs this file (`download_genome`) before running `uncompress_genome`!

Try: `snakemake -p uncompress_genome`

You'll see that it runs two rules: first the `download_genome` rule, then the `uncompress_genome` rule.
![](../../images/snakemake_twosteps.jpeg)

!!! recap
    - `input` and `output` (and other things) can be in any order, as long as they are before `shell`
    - for each of the above elements, their contents can be all on one line, or form a block by indenting
    - you can make lists for multiple input or output files by separating filenames with a comma
    - rule names can be any valid variable, which basically means letters and underscores; you can use numbers after a first character; no spaces!

### Continue decorating

The remaining rules need inputs and outputs so we can link up all the rules thru to the final variant calling step.

!!! tip
    The [second vidlet](https://video.ucdavis.edu/media/snakemake+intro+2+try+1/0_t1dpuzly) and [third vidlet](https://video.ucdavis.edu/media/snakemake+intro+3+try+1/0_gwnss4kq) cover the remaining content by adding more detail to the Snakefile and wrapping up with the final workflow rule.

Run the rules one at a time to figure out what the output files are. Check the file time stamps (`ls -lht`) to track more recent output files. The inputs are specified in the shell command section.

- Note that you can do a dry run of the rule with `snakemake -n <rule name>` to check how snakemake is interpreting the rule input(s), output(s), and shell command(s), without actually running the command(s) or creating any output(s).

Sometimes a command requires multiple input files but only explicitly states one (the software assumes that if a certain file exists the other required files must exist). While the command doesn't have a space for it, we will define all inputs in the rule's `input:` section.

Try filling in `input` and `output` sections. Watch the vidlets if you are struggling with this!

A complete Snakefile is [here](./example_snakefile.md). Note that there are many ways to concisely enter the input and output files and this is just one example version!

### Running lots of rules all at once
Once you've fixed the rules `index_genome_bwa` and `map_reads`, you should be able to run everything up to the rule `index_genome_samtools` by running:
```
(snaketest) $ snakemake -p index_genome_samtools
```

This also serves as a good way to check that you have all the correct input/output information. You'll have files left over if you forgot to put them in output.

### Re-running everything
You can tell snakemake to delete all the inputs/outputs for a particular rule (including preceding rules) by running:
```
(snaketest) $ snakemake --delete-all-output index_genome_samtools
```

### Using filenames instead of rule names
You don't actually need to use the rule names (this will be important later on!). You can tell snakemake what file you want produced and it will run all the rules necessary to get to that output (as long as the upstream rules have been properly linked!).

```
(snaketest) $ snakemake -p SRR2584857.sam
```

will also work to run the rule `map_reads`, but you don't have to remember the rule name (which can be arbitrary).

### Some python shortcuts
There are several ways to more efficiently and cleanly specify inputs/outputs in the Snakefile. Here are some shortcuts:

Take the `uncompress_genome` rule we decorated above:
```
rule uncompress_genome:
    input: "ecoli-rel606.fa.gz"
    shell:
        "gunzip ecoli-rel606.fa.gz"
```
It can also be written as follows with wildcards `{input}` and `{output}`. Wildcards operate entirely within a single rule, not across rules. This means that we can use different definitions for `{input}` and `{output}` for each rule and they won't conflict.
```
rule uncompress_genome:
    input: "ecoli-rel606.fa.gz"
    output: "ecoli-rel606.fa"
    shell:
        "gunzip -c {input} > {output}"
```

Multiple inputs files can be separated by commas and written on their own lines. The input files can be assigned variable names that are accessed in the `shell:` block with `input.<input file variable>`. The `\` tells the shell that this is one command written over two lines in the file. This rule also demonstrates an input file (.bai) that is not explicitly stated but need for the shell command.
```
rule samtools_mpileup:
    input:
        index="ecoli-rel606.fa",
        sorted="SRR2584857.sorted.bam",
        sorted_bai="SRR2584857.sorted.bam.bai"
    output:
        "variants.raw.bcf"
    shell:
        """
        samtools mpileup -u -t DP -f {input.index} {input.sorted} | \
        bcftools call -mv -Ob -o - > {output}
        """
```

Since the Snakefile is written in Python, we can use Python functions! Here, we are consolidating the list of .fa files into a single line of code. The expansion (`exp`) tells Python that there is a common file name pattern (`ecoli-rel606.fa.`) with different endings (`{ext}`) that are specified using a list (`ext=['sa', 'amb', 'ann', 'pac', 'bwt']`).
```
rule map_reads:
    input:
        genome = "ecoli-rel606.fa",
        reads = "SRR2584857_1.fastq.gz",
        idxfile = expand("ecoli-rel606.fa.{ext}", ext=['sa', 'amb', 'ann', 'pac', 'bwt'])
    output: "SRR2584857.sam"
    shell:
        "bwa mem -t 4 {input.genome} {input.reads} > {output}"
```

## Run the whole pipeline!

Once you've decorated the entire Snakefile, you should be able to run through the whole workflow:
```
# run final rule of workflow
(snaketest) $ snakemake -p make_vcf

# or run final output of final rule
(snaketest) $ snakemake -p variants.vcf
```

Alternatively, you can create a default rule called `all` at the top of the Snakefile:
```
rule all:
    input: "variants.vcf"
```
and run:
```
(snaketest) $ snakemake -p all
```

This rule runs through the entire workflow with a single command! This is much better than running each command one by one!

!!! recap
    - Snakefile defines a snakemake workflow
    - the rules specify steps in the workflow
    - at the moment (and in general), they run shell commands
    - you can "decorate" the rules to tell snakemake how they depend on each other
    - this decoration comes in the form of `input:` and `output:` lists of one or more files, quoted, separated by commas
    - this is how you connect rules: by saying which rules take which files as inputs and/or produce what outputs
    - snakemake cares about tabs

### Looking at VCF files
Finally, let's look at the output!

```
(snaketest) $ less variants.vcf
```

We can look at the alignment by running the [alignment viewer](http://samtools.sourceforge.net/tview.shtml) in samtools:
```
(snaketest) $ samtools tview -p ecoli:4202391 SRR2584857.sorted.bam ecoli-rel606.fa
```

!!! tip "Navigating samtools tview"
    - Access the help menu in `samtools tview` by hitting `shift` key and `?` key.
    - Exit help menu by hitting `q` key
    - Go to a specific position location in the alignment by hitting `/` key.

### Wrapping up
Exit conda environment to return to base environment:
```
(snaketest) $ conda deactivate
```

Hopefully, you have now:

!!! recap

    - learned how to write basic workflows with snakemake rules
    - learned variable substitution for snakemake rules
    - learned wildcard matching for snakemake rules
    - understood why workflow systems can help you do your computing more easily!
