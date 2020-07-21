# download data from https://osf.io/vzfc6/

rule download_data:
    shell:
        "wget https://osf.io/4rdza/download -O SRR2584857_1.fastq.gz"

rule download_genome:
    shell:
        "wget https://osf.io/8sm92/download -O ecoli-rel606.fa.gz"

rule uncompress_genome:
    shell:
        "gunzip ecoli-rel606.fa.gz"

rule index_genome_bwa:
    shell:
        "bwa index ecoli-rel606.fa"

rule map_reads:
    shell:
        "bwa mem -t 4 ecoli-rel606.fa SRR2584857_1.fastq.gz > SRR2584857.sam"

rule index_genome_samtools:
    shell:
        "samtools faidx ecoli-rel606.fa"

rule samtools_import:
    shell:
        "samtools import ecoli-rel606.fa.fai SRR2584857.sam SRR2584857.bam"

rule samtools_sort:
    shell:
        "samtools sort SRR2584857.bam -o SRR2584857.sorted.bam"

rule samtools_index_sorted:
    shell: "samtools index SRR2584857.sorted.bam"


rule samtools_mpileup:
    shell:
        """samtools mpileup -u -t DP -f ecoli-rel606.fa SRR2584857.sorted.bam | \
    bcftools call -mv -Ob -o - > variants.raw.bcf"""

rule make_vcf:
    shell: "bcftools view variants.raw.bcf > variants.vcf"

# at end, run:
## samtools tview -p ecoli:4202391 SRR2584857.sorted.bam ecoli-rel606.fa
