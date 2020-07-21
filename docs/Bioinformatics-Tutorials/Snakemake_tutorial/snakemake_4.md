# Continue Decorating

The remaining rules need inputs and outputs so we can link up all the rules through to the final variant calling step. The Part 2 and Part 3 of video tutorials cover the remaining content by adding more detail to the Snakefile and wrapping up with the final workflow rule.

Part 2: Continue decorating the Snakefile

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=0_t1dpuzly&flashvars[mediaProtocol]=rtmp&amp;flashvars[streamerType]=rtmp&amp;flashvars[streamerUrl]=rtmp://www.kaltura.com:1935&amp;flashvars[rtmpFlavors]=1&amp;flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=0_01u9hxvk" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>

Part 3: Running through the entire Snakemake workflow

<iframe id="kaltura_player" src="https://cdnapisec.kaltura.com/p/1770401/sp/177040100/embedIframeJs/uiconf_id/29032722/partner_id/1770401?iframeembed=true&playerId=kaltura_player&entry_id=0_gwnss4kq&flashvars[mediaProtocol]=rtmp&amp;flashvars[streamerType]=rtmp&amp;flashvars[streamerUrl]=rtmp://www.kaltura.com:1935&amp;flashvars[rtmpFlavors]=1&amp;flashvars[localizationCode]=en&amp;flashvars[leadWithHTML5]=true&amp;flashvars[sideBarContainer.plugin]=true&amp;flashvars[sideBarContainer.position]=left&amp;flashvars[sideBarContainer.clickToClose]=true&amp;flashvars[chapters.plugin]=true&amp;flashvars[chapters.layout]=vertical&amp;flashvars[chapters.thumbnailRotator]=false&amp;flashvars[streamSelector.plugin]=true&amp;flashvars[EmbedPlayer.SpinnerTarget]=videoHolder&amp;flashvars[dualScreen.plugin]=true&amp;flashvars[Kaltura.addCrossoriginToIframe]=true&amp;&wid=0_kjfuqewn" width="608" height="402" allowfullscreen webkitallowfullscreen mozAllowFullScreen allow="autoplay *; fullscreen *; encrypted-media *" sandbox="allow-forms allow-same-origin allow-scripts allow-top-navigation allow-pointer-lock allow-popups allow-modals allow-orientation-lock allow-popups-to-escape-sandbox allow-presentation allow-top-navigation-by-user-activation" frameborder="0" title="Kaltura Player"></iframe>

Run the rules one at a time to figure out what the output files are. Check the file time stamps (`ls -lht`) to track more recent output files. The inputs are specified in the shell command section.

!!! Tip
    You can do a dry run of the rule with `snakemake -n <rule name>` to check how Snakemake is interpreting the rule input(s), output(s), and shell command(s), without actually running the command(s) or creating any output(s).

Sometimes a command may require multiple input files but only explicitly state one in the command (the software assumes that if a certain file exists the other required files must exist). To avoid potential errors in rule dependencies, we will define all inputs in the rule's `input:` section.

In this workflow, the rule `map_reads` is a good example of such a behavior. `bwa` is used to generate the mapped reads (`.sam`) file and requires the reference genome and the index files without explicitly referring to the index files in the command. We can add an additional input variable for index files in `map_reads` to define all input files to Snakemake:

```
# Map the raw reads to the reference genome
rule map_reads:
    input:
        genome = "ecoli-rel606.fa",
        reads = "SRR2584857_1.fastq.gz",
        idxfile = expand("ecoli-rel606.fa.{ext}", ext=['sa', 'amb', 'ann', 'pac', 'bwt'])
    output: "SRR2584857.sam"
    shell:
        "bwa mem -t 4 {input.genome} {input.reads} > {output}"
```        

Follow along the video tutorials to fill in the `input` and `output` sections for the remaining rules.

The complete Snakefile is [here](./example_snakefile.md). Note that there are many ways to concisely enter the input and output files and this is just one example version!

### Running lots of rules all at once

Once you've fixed the rules `index_genome_bwa` and `map_reads`, you should be able to run everything up to the rule `index_genome_samtools` by running:

```
(snaketest) $ snakemake -p index_genome_samtools
```

This also serves as a good way to check that you have all the correct input/output information. You'll have files left over if you forgot to put them in output.

### Re-running rules

Snakemake also has the option to delete all the inputs/outputs for a particular rule (including preceding rules), shown here by running the `index_genome_samtools` rule:
```
(snaketest) $ snakemake --delete-all-output index_genome_samtools
```

### Using filenames instead of rule names

You don't actually need to use the rule names *(this will be important later on!)*. Instead of rule names, you can specify the required output file in Snakemake which will trigger execution of all the upstream linked rules necessary to produce the file.

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

Multiple inputs files can be separated by commas and written on their own lines. The input files can be assigned variable names that are accessed in the `shell:` block with `input.<input file variable>`. The `\` tells the shell that this is one command written over two lines in the file. Also, similar to `map_reads`, the `samtools_mpileup` rule also includes an input file (`.bai`) that is not explicitly stated but required to complete the command.

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

Since the Snakefile is written in Python, we can also use Python functions! As an example, we will consolidate the list of reference genome index files into a single line of code. The expansion (`expand`) tells Python that there is a common file name pattern (`ecoli-rel606.fa.`) with different endings (`{ext}`) that are specified using a list (`ext=['sa', 'amb', 'ann', 'pac', 'bwt']`).

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

Once you've decorated the entire Snakefile, you should be able to run through the whole workflow using either the rule name:

```
(snaketest) $ snakemake -p make_vcf
```

or by using the file name:

```
(snaketest) $ snakemake -p variants.vcf
```

Recollect that Snakemake will execute the first rule in the Snakefile by default. We can use this feature by creating a default rule called `all` at the top of the Snakefile:

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
    - Snakefile defines a Snakemake workflow
    - the rules specify steps in the workflow
    - at the moment (and in general), they run shell commands
    - you can "decorate" the rules to link the dependencies between rules
    - basic decoration comes in the form of `input:` and `output:` which are lists of one or more files, quoted, separated by commas
    - the rules are connected by matching filenames
    - tabs are important syntax feature in Snakemake

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

- learned how to write basic workflows with Snakemake rules
- learned variable substitution for Snakemake rules
- learned wildcard matching for Snakemake rules
- understood why workflow systems can help you do your computing more easily!
