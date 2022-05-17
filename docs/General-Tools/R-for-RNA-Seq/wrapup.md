---
title: Wrap-up
---

Files for this workshop can be downloaded from GitHub at https://github.com/nih-cfde/training-rstudio-binder.
Specifically, scripts and data are in the `GTEx` directory.

### References

-   [R for Data Science by Hadley Wickham and Garrett
    Grolemund](https://r4ds.had.co.nz/index.html)
-   [Rouillard et al. 2016. The Harmonizome: a collection of processed
    datasets gathered to serve and mine knowledge about genes and
    proteins. Database
    (Oxford).](http://database.oxfordjournals.org/content/2016/baw100.short)
-   [GTEx Data Portal](https://gtexportal.org/home/)

### Additional Resources

-   [RStudio cheat sheet for
    readr](https://raw.githubusercontent.com/rstudio/cheatsheets/master/data-import.pdf)
-   [RStudio cheat sheet for
    dplyr](https://raw.githubusercontent.com/rstudio/cheatsheets/master/data-transformation.pdf)
-   [RStudio cheat sheet for data Wrangling with
    dplyr](https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf)
-   [ggplot point
    shapes](http://www.sthda.com/english/wiki/ggplot2-point-shapes)
-   [Angus 2019 Intro to R
    Lesson](https://angus.readthedocs.io/en/2019/R_Intro_Lesson.html)
-   [Angus 2019 Differential Gene Expression in R
    Lesson](https://angus.readthedocs.io/en/2019/diff-ex-and-viz.html)
-   [Software Carpentry R
    Lesson](http://swcarpentry.github.io/r-novice-inflammation/)

### Appendix

Below is a consolidated list of the commands used in this workshop. 

```r
## Introduction

2 + 2 * 100
log10(0.05)

pval <- 0.05
pval

-log10(pval)


favorite_genes <- c("BRCA1", "JUN",  "GNRH1", "TH", "AR")
favorite_genes

#install.packages("ggplot2")

library(ggplot2)
library(tidyr)
library(dplyr)


## Importing Data

samples <- read.csv("./data/samples.csv")

#View(samples)
head(samples)
tail(samples)
str(samples)
summary(samples)


counts <- read.csv("./data/countData.HEART.csv", row.names = 1)
dim(counts)
head(counts)[1:5]

results <- read.table("./data/GTEx_Heart_20-29_vs_50-59.tsv")
head(results)

dim(samples)

dplyr::count(samples, SMTS) 

head(dplyr::count(samples, SMTS, SEX))

head(dplyr::count(samples, SMTS, SEX, AGE, DTHHRDY ) )

## Visualizing Data

ggplot(samples, aes(x = SMTS)) +
  geom_bar(stat = "count")


ggplot(samples, aes(x = SMTS)) +
  geom_bar(stat = "count") + 
  coord_flip()


ggplot(samples, aes(x = SMTS, color = AGE)) +
  geom_bar(stat = "count") + 
  coord_flip()


ggplot(samples, aes(x = SMTS, fill = AGE)) +
  geom_bar(stat = "count") + 
  coord_flip()


ggplot(samples, aes(x = SMTS, fill = AGE)) +
  geom_bar(stat = "count") + 
  coord_flip() +
  facet_wrap(~SEX)


ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point() 


ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05))


ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ifelse( adj.P.Val < 0.05, "p < 0.05", "NS"))) +
  geom_hline(yintercept = -log10(0.05)) 


ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ifelse( adj.P.Val < 0.05, "p < 0.05", "NS"))) +
  geom_hline(yintercept = -log10(0.05))  +
  theme(legend.position = "bottom") +
  labs(color = "20-29 vs 50-59 year olds", 
       subtitle = "Heart Tissue Gene Expression")


ggplot(samples, aes(x = SMCENTER, y = SMRIN)) +
  geom_boxplot() +
  geom_jitter(aes(color = SMRIN))


## Wrangle Data

results %>% 
  filter(adj.P.Val < 0.05) %>% 
  head()


results %>% 
  filter(logFC > 1 | logFC < -1) %>%
  head()


results %>% filter(adj.P.Val < 0.05,
                   logFC > 1 | logFC < -1) %>%
  arrange(adj.P.Val) %>%
  head()



resultsDEGs <- results %>% filter(adj.P.Val < 0.05,
                                  logFC > 1 | logFC < -1) %>%
  arrange(adj.P.Val) %>% 
  rownames(.)
resultsDEGs


colData <- read.csv("./data/colData.HEART.csv", row.names = 1)
head(colData)

head(rownames(colData) == colnames(counts))
head(colnames(counts))
head(rownames(colData))


colData_tidy <-  colData %>%
  mutate(SAMPID = gsub("-", ".", SAMPID))  
rownames(colData_tidy) <- colData_tidy$SAMPID

mycols <- rownames(colData_tidy)
head(mycols)


counts_tidy <- counts %>%
  select(all_of(mycols))

head(rownames(colData_tidy) == colnames(counts_tidy))


genes <- read.table("./data/ensembl_genes.tsv", sep = "\t",  header = T, fill = T)
head(genes)


resultsSymbol <- results %>%
  mutate(name = row.names(.))
head(resultsSymbol)


resultsName <- left_join(resultsSymbol, genes, by = "name")
head(resultsName)


resultsNameTidy <- resultsName %>%
  filter(adj.P.Val < 0.05,
         logFC > 1 | logFC < -1) %>%
  arrange(adj.P.Val) %>%
  select(name, description, id, logFC, AveExpr, adj.P.Val)
head(resultsNameTidy)


resultsNameTidyIds <- resultsNameTidy %>%
  drop_na(id) %>%
  pull(id)
resultsNameTidyIds


counts_tidy_slim <- counts_tidy %>%
  mutate(id = row.names(.)) %>%
  filter(id %in% resultsNameTidyIds)
dim(counts_tidy_slim)
head(counts_tidy_slim)[1:5]
tail(counts_tidy_slim)[1:5]

counts_tidy_long <- counts_tidy_slim %>%
  pivot_longer(cols = all_of(mycols), names_to = "SAMPID", 
               values_to = "counts") 
head(counts_tidy_long)


counts_tidy_long_joined <- counts_tidy_long%>%
  inner_join(., colData_tidy, by = "SAMPID") %>%
  inner_join(., genes, by = "id") %>%
  arrange(desc(counts))
head(counts_tidy_long_joined)

library(scales)

counts_tidy_long_joined %>%
  ggplot(aes(x = AGE, y = counts)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~name, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust  = 1),
        strip.text = element_text(face = "italic")) +
  scale_y_log10(labels = label_number_si()) 
```
