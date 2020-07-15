# Example Snakefile

    #Data downloaded from: https://osf.io/vzfc6/

    # Default rule to run entire workflow, only works once inputs/outputs correctly filled in all rules
    rule all:
        input: "variants.vcf"

    # Download raw sequence data
    rule download_data:
        output: "SRR2584857_1.fastq.gz"
        shell:
            "wget https://osf.io/4rdza/download -O {output}"

    # Download reference genome
    rule download_genome:
        output: "ecoli-rel606.fa.gz"
        shell:
            "wget https://osf.io/8sm92/download -O {output}"

    # The .gz file must be uncompressed before it can be accessed
    rule uncompress_genome:
        input: "ecoli-rel606.fa.gz"
        output: "ecoli-rel606.fa"
        shell:
            "gunzip -c {input} > {output}"

    # Create an index for the reference genome for bwa
    rule index_genome_bwa:
        input: "ecoli-rel606.fa"
        output:
            expand("ecoli-rel606.fa.{ext}", ext=['sa', 'amb', 'ann', 'pac', 'bwt'])
        shell:
            "bwa index {input}"

    # Map the raw reads to the reference genome
    rule map_reads:
        input:
            genome = "ecoli-rel606.fa",
            reads = "SRR2584857_1.fastq.gz",
            idxfile = expand("ecoli-rel606.fa.{ext}", ext=['sa', 'amb', 'ann', 'pac', 'bwt'])
        output: "SRR2584857.sam"
        shell:
            "bwa mem -t 4 {input.genome} {input.reads} > {output}"

    # Create an index for the reference genome fasta file for samtools
    # Note: this indexing step is different and unrelated to the one above for bwa
    rule index_genome_samtools:
        input: "ecoli-rel606.fa"
        output: "ecoli-rel606.fa.fai"
        shell:
            "samtools faidx {input}"

    # Convert .sam to .bam file
    rule samtools_import:
        input:
            index="ecoli-rel606.fa.fai",
            samfile="SRR2584857.sam"
        output:"SRR2584857.bam"
        shell:
            """
            samtools view -bt {input.index} -o {output} {input.samfile}
            """
            # original command with samtools v1.9
            ## samtools import {input.index} {input.samfile} {output}
            # but it gave segmentation fault error with samtools v1.10, so here we use samtools view instead

    # Sort the bam alignment file
    rule samtools_sort:
        input: "SRR2584857.bam"
        output: "SRR2584857.sorted.bam"
        shell:
            "samtools sort {input} -o {output}"

    # Create an index for the sorted bam file
    # again, this is different from the above indexing steps
    rule samtools_index_sorted:
        input: "SRR2584857.sorted.bam"
        output: "SRR2584857.sorted.bam.bai"
        shell: "samtools index {input}"

    # Generate pileup file with samtools, then call variants with bcftools
    # From samtools doc: 'Pileup format consists of TAB-separated lines, with each line representing the pileup of reads at a single genomic position.'
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

    # Convert bcf to vcf file so we can look at the resulting mappings
    rule make_vcf:
        input: "variants.raw.bcf"
        output: "variants.vcf"
        shell: "bcftools view {input} > {output}"

    # at end, run this command in the terminal to view the mapping:
    ## samtools tview -p ecoli:4202391 SRR2584857.sorted.bam ecoli-rel606.fa
    
