This is the WDL script for our demo workflow. It should be saved with the ".wdl" extension.

```
workflow GWAS {
    String sample_name
    File inputvcf
    File inputpheno

    call run_vcftools {
      input:
        sample_name = sample_name,
        inputvcf = inputvcf
    }
    call plink_missing_rates {
      input:
        inputpheno = inputpheno,
        inputminoralleles = run_vcftools.outminoralleles,
        inped = run_vcftools.outped,
        inmap = run_vcftools.outmap
    }
    call plink_binary {
      input:
        sample_name = sample_name,
        inped = run_vcftools.outped,
        inmap = run_vcftools.outmap
    }
    call plink_association {
      input:
        sample_name = sample_name,
        inputpheno = inputpheno,
        inputminoralleles = run_vcftools.outminoralleles,
        inputbed = plink_binary.outbed,
        inputfam = plink_binary.outfam,
        inputbim = plink_binary.outbim
    }
    call run_R {
      input:
        sample_name = sample_name,
        inadj = plink_association.outadj,
        inassoc = plink_association.outassoc
     }
}

task run_vcftools {
    String sample_name
    File inputvcf

    command <<<
      vcftools --vcf ${inputvcf} --plink --out ${sample_name}
      cat ${inputvcf} | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{{if($3==".")$3=$1":"$2;}print $3,$5;}' > minor_alleles
    >>>

    output {
      File outmap = '${sample_name}.map'
      File outped = '${sample_name}.ped'
      File outminoralleles = 'minor_alleles'
    }
    # the coatColor.log file is going to stderr; do not define as an output - get an error because it can't find the output file name

    runtime {
      docker: 'mlim13/demo_gwas:tag0'
    }
}

task plink_missing_rates {
    File inputpheno
    File inputminoralleles
    File inped
    File inmap

    command <<<
      ## must state input files explicitly (can use --file OR --ped and --map, not both flags)
      ## cannot use '\' for line breaks it seems
      plink --ped ${inped} --map ${inmap} --make-pheno ${inputpheno} "yellow" --missing --out miss_stat --noweb --dog --reference-allele ${inputminoralleles} --allow-no-sex --adjust
    >>>

    output {
      File outimiss = 'miss_stat.imiss'
      File outlmiss = 'miss_stat.lmiss'
      File outstatno = 'miss_stat.nosex'
    }

    runtime {
      ## this is a docker specifically for plink.
      docker: 'gelog/plink:latest'
    }
}

task plink_binary {
    String sample_name
    File inped
    File inmap

    command {
      plink --ped ${inped} --map ${inmap} --allow-no-sex --dog --make-bed --noweb --out ${sample_name}.binary
    }

    output {
      File outfam = '${sample_name}.binary.fam'
      File outbed = '${sample_name}.binary.bed'
      File outbim = '${sample_name}.binary.bim'
      File outno = '${sample_name}.binary.nosex'
    }

    runtime {
      ## this is a docker specifically for plink.
      docker: 'gelog/plink:latest'
    }
}

task plink_association {
    String sample_name
    File inputminoralleles
    File inputpheno
    File inputfam
    File inputbed
    File inputbim

    command <<<
      plink --bed ${inputbed} --bim ${inputbim} --fam ${inputfam} --make-pheno ${inputpheno} "yellow" --assoc --reference-allele ${inputminoralleles} --allow-no-sex --adjust --dog --noweb --out ${sample_name}
    >>>

    output {
      File outadj = '${sample_name}.assoc.adjusted'
      File outassoc = '${sample_name}.assoc'
      File outno = '${sample_name}.nosex'
    }

    runtime {
      ## this is a docker specifically for plink.
      docker: 'gelog/plink:latest'
    }
}

task run_R {
    String sample_name
    File inadj
    File inassoc

    command <<<
      unad_cutoff_sug=$(tail -n+2 ${inadj} | awk '$10>=0.05' | head -n1 | awk '{print $3}')
      unad_cutoff_conf=$(tail -n+2 ${inadj} | awk '$10>=0.01' | head -n1 | awk '{print $3}')
      Rscript -e 'args=(commandArgs(TRUE));library(qqman);'\ 'data=read.table("${inassoc}", header=TRUE); data=data[!is.na(data$P),];'\ 'bitmap("${sample_name}_man.bmp", width=20, height=10);'\ 'manhattan(data, p = "P", col = c("blue4", "orange3"),'\ 'suggestiveline = 12,'\ 'genomewideline = 15,'\ 'chrlabs = c(1:38, "X"), annotateTop=TRUE, cex = 1.2);'\ 'graphics.off();' $unad_cutoff_sug $unad_cutoff_conf
    >>>

    output {
      File outbmp = '${sample_name}_man.bmp'
    }

    runtime {
      docker: 'mlim13/demo_gwas:tag0'
    }
}
```
