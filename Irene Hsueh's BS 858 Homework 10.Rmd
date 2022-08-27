---
title: "Irene Hsueh's BS 858 Homework 10"
author: "Irene Hsueh"
date: "11/29/2021"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
```

### 
```{r}
#Reading in CSV File
snp_info <- read.csv("C:/Irene Hsueh's Documents/MS Applied Biostatistics/BS 858 - Statistical Genetics I/Class 11 - Association Analysis of Rare Variants/Homework 10/snp_info_hw10.csv") %>%
#Selecting and Renaming Variables  
  dplyr::select(variant_name = varname, 
                gene, 
                mutation = vartype, 
                variant)


#Reading in CSV File
genotypes <- read.csv("C:/Irene Hsueh's Documents/MS Applied Biostatistics/BS 858 - Statistical Genetics I/Class 11 - Association Analysis of Rare Variants/Homework 10/pheno_geno_hw10.csv") %>%
#Selecting and Renaming Variables  
  dplyr::select(var1:var21, 
                age, sex, 
                population1 = pop1, 
                population2 = pop2, 
                population3 = pop3, 
                quantitative_trait = qt, 
                affected)
head(snp_info, 10)
head(genotypes, 10)
```



```{r}
#Minor Allele Frequencies 
snp_info$maf <- apply(genotypes[1:21], 2, function(x) sum(x)/(2*length(x)))


#Rare SNPs
snp_info$rare_variants <- ifelse(snp_info$maf <= 0.01, yes=TRUE, no=FALSE)
rare_index <- which(snp_info$rare_variants==TRUE)
common_index <- which(snp_info$rare_variants==FALSE)
genotypes_rare <- genotypes[, -common_index]
```




### Weighted Anaysis
```{r}
weights <- 1/sqrt(nrow(genotypes) * snp_info$maf * (1-snp_info$maf) )

#Weights for All Variants
weights_matrix_all <- sapply(1:21, function(x) genotypes[,x]*weights[x])

genotypes$weighted_score_all <- apply(weights_matrix_all, 1, sum)
summary(genotypes$weighted_score_all)
sd(genotypes$weighted_score_all)


#Weights for Rare Variants
weights_matrix_rare <- sapply(1:length(rare_index), function(x) genotypes_rare[,x]*weights[x])

genotypes_rare$weighted_score_rare <- apply(weights_matrix_rare, 1, sum)
summary(genotypes_rare$weighted_score_rare)
sd(genotypes_rare$weighted_score_rare)


#Weighted Analysis using All Variants
weighted_model_all_variants <- lm(quantitative_trait ~ weighted_score_all + sex + age + population1 + population2 + population3, data=genotypes)
summary(weighted_model_all_variants)

#Weighted Analysis using Rare Variants
weighted_model_rare_variants <- lm(quantitative_trait ~ weighted_score_rare + sex + age + population1 + population2 + population3, data=genotypes_rare)
summary(weighted_model_rare_variants)
```



### Burden Tests
```{r}
#Nonsynonymous SNPs
nonysnonymous_index <- which(snp_info$mutation=="Nonsynonymous")
other_mutations_index <- which(snp_info$mutation!="Nonsynonymous")
genotypes_nonsynonymous <- genotypes[, -other_mutations_index]


#Nonsynonymous & Rare SNPs
nonsynonymous_rare_index <- which(snp_info$rare_variants==TRUE & snp_info$mutation=="Nonsynonymous")
other_common_mutations_index <- which(snp_info$rare_variants==FALSE | snp_info$mutation!="Nonsynonymous")
genotypes_rare_nonsynonymous <- genotypes[, -other_mutations_index]


#CMC Burden Test
genotypes_rare$cmc <- apply(genotypes_rare[1:length(rare_index)], 1, sum)
table(genotypes_rare$cmc)

genotypes_nonsynonymous$cmc <- apply(genotypes_nonsynonymous[1:length(nonysnonymous_index)], 1, sum)
table(genotypes_nonsynonymous$cmc)

genotypes_rare_nonsynonymous$cmc <- apply(genotypes_rare_nonsynonymous[1:length(nonsynonymous_rare_index)], 1, sum)
table(genotypes_rare_nonsynonymous$cmc)



#Weighted Analysis using Rare Variants
cmc_model_rare <- lm(quantitative_trait ~ cmc + sex + age + population1 + population2 + population3, data=genotypes_rare)
summary(cmc_model_rare)

#Weighted Analysis using Nonsynonymous Variants
cmc_model_nonsynonymous <- lm(quantitative_trait ~ cmc + sex + age + population1 + population2 + population3, data=genotypes_nonsynonymous)
summary(cmc_model_nonsynonymous)

#Weighted Analysis using Rare & Nonsynonymous Variants
cmc_model_nonsynonymous_rare <- lm(quantitative_trait ~ cmc + sex + age + population1 + population2 + population3, data=genotypes_rare_nonsynonymous)
summary(cmc_model_nonsynonymous_rare)

```


### More CMC Tests
```{r}
#Weighted Analysis using Rare & Nonsynonymous Variants
cmc_model_nonsynonymous_rare2 <- lm(quantitative_trait ~ cmc + sex + age, data=genotypes_rare_nonsynonymous)
summary(cmc_model_nonsynonymous_rare2)
```


