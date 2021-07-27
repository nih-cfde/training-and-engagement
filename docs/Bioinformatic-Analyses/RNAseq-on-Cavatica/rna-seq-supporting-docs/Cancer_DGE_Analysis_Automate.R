~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Author: Saranya Canchi
  # Filename: Cancer_DGE_Analysis_Automate.R
  # Purpose: R script for differential gene expression between 
  #          pediatric cancer types using DESeq2 
  # Version: 1.0
  # Date: 01/22/2021
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
# Install libraries####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("tximport",
                       "regionReport",
                       "BiocStyle"), 
                       update=FALSE,
                       ask=FALSE)

# Load libraries####

library(GenomicFeatures)
library(DESeq2)
library(tximport)
library(regionReport)
library(ggplot2)
library(knitr)

# Generate transcript and gene name table for gene level summary####
## Using GenomicFeatures pkg to read in the reference GTF files - following steps from DESeq2 app

data_dir <- "/sbgenomics/project-files"
txdb <- makeTxDbFromGFF(file= file.path(data_dir,"gencode.v27.annotation.gtf"))
k <- keys(txdb, keytype = "TXNAME")

## HGNC gene name not available in list of filters. Using Ensemble gene name for one to one mapping. 

tx2gene <- select(txdb, k, "GENEID", "TXNAME") 

# Read the phenotype data####

pheno_data <- read.csv(file.path(data_dir,"phenotype_filtered.csv"))

## Converting covariates of interest to factors
### Setting Ependymoma as reference factor level for histology variable

pheno_data$histology <- factor(pheno_data$histology, levels=c("Ependymoma", "Medulloblastoma"))
pheno_data$tumor_location <- factor(pheno_data$tumor_location)
pheno_data$diagnosis_age_range <- factor(pheno_data$diagnosis_age_range)

# Generate gene level summary using tximport####

head(list.files(data_dir))
files <- file.path(data_dir, pheno_data$name)
names(files) <- pheno_data$sample_id
txi_sum <- tximport(files,
                    type="kallisto",
                    tx2gene=tx2gene,
                    ignoreTxVersion = TRUE)

# DESeq2 import and analysis####

dds_cancer <- DESeqDataSetFromTximport(txi_sum,
                                       colData = pheno_data,
                                       design = ~ diagnosis_age_range + tumor_location + histology)

# DGE analysis#### 

dds_cancer <- DESeq(dds_cancer, 
                    fit="parametric", 
                    test='Wald', 
                    betaPrior=TRUE, 
                    parallel=TRUE)

res_cancer <- results(dds_cancer,
                      contrast=c("histology","Medulloblastoma","Ependymoma"),
                      alpha=0.05)

## Write the results to a table

output_dir <- "/sbgenomics/output-files"
res_cancer_order <- res_cancer[order(res_cancer$pvalue),]
write.csv(as.data.frame(res_cancer_order),
          file=file.path(output_dir,"Cancer_DESeq2_DGE_results.csv"))

## Write normalized counts to a file

norm_counts <- counts(dds_cancer,normalized=TRUE)
write.table(norm_counts, 
            file=file.path(output_dir,"Cancer_DESeq2_normalized_counts.txt"),
            sep="\t",
            quote=FALSE,
            col.names = NA) 

# Create a report with all the visualization from the DESeq2 vignette####

dir.create(file.path(output_dir,"DESeq2-Report"), 
           showWarnings = FALSE, 
           recursive = TRUE)

report <- DESeq2Report(dds = dds_cancer, 
                       project = 'Pediatric Cancer DGE Analysis with DESeq2', 
                       intgroup = c('histology','diagnosis_age_range'),
                       res = res_cancer,   
                       outdir = file.path(output_dir,"DESeq2-Report"),
                       output = 'Cancer-DESeq2-Report', 
                       theme = theme_bw())

save.image(file=file.path(output_dir,paste0("Cancer_DGE_",Sys.Date(), ".env.Rdata")))

