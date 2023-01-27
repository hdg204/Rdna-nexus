# Rdna-nexus

## Introduction

This repository contains a few useful functions for analysing genetic data from the UK Biobank using the RStudio Workbench implementation on DNA Nexus. These functions were written for use on the University of Exeter's project, but as none of our curated files are used, this should run universally. There are currently three functions included:

* **extract_snp**, which is for extracting everbody's genotype for a particular SNP
* **extract_snp_bulk**, which is for extracting everbody's genotype for lots of SNPs
* **generate_grs**, which combines lots of SNPs into a Genetic Risk Score

## Installation

To make all functions available, as well as an example file, use
```
library(devtools)
source_url("https://raw.githubusercontent.com/hdg204/Rdna-nexus/main/install.R")
```

This will source all R functions in the repository.

## extract_snp

*Extract the data of one SNP from the imputed genotype file*

This function takes the chromosome and base pair position of one SNP and extracts the genotype data. The inputs are two numeric variables for base pair and position.

E.g. to extract the data for the SNP on chromosome 8, position 38644914, run:

`a=extract_snp(8,38644914)`

The output of this function is a list with two elements.

`a$snp_data` is a dataframe with two headings, one for the eid, and one for the genotype of everyone in Biobank

`a$variant_info` contains information about the SNP itself: a dataframe with headings chromosome, position, rsid, number_of_alleles, allele0, allele1
where the rsid is extracted from the bgen file.

![image](https://user-images.githubusercontent.com/36624710/215060277-b734c84f-5708-4a3b-aa52-82957eb531c0.png)

An image of the genotype data cannot be provided as this contains individual-level genetic data.

## extract_snp_bulk

*Extract the data of many SNPs*

While many SNPs can be extracted using `extract_snp`, this is slow as it would need to read the bgen for each SNP. The most efficient way is to read the bgen for each chromosome, since these are split in the UKBB data. This function is built for extracting multiple SNPs from a file.

The file should be a tab separated file, with the first two columns being chromosome and base pair position. An example of this is given in https://github.com/hdg204/Rdna-nexus/blob/main/Example_GRS, where only the first two columns are used.

After installation, this function can be run using:

`a=extract_snp_bulk('Example_GRS')`

It can be a little slow but it prints to screen which chromosome it's on. The output is similar to extract_snp. It is a list of two elements, genotypes and variants. This time there are multiple columns for the genotypes, named 'chr:bp_A0_A1', and multiple rows for the variants, one for each SNP. The variants table looks like hte following:

![image](https://user-images.githubusercontent.com/36624710/215066013-9689fd7f-5bae-447e-b16b-407d88a91397.png)


## generate_grs

*Calculate a genetic risk score using SNP weights*

This function takes a weight file, which must have the columns chromosome, base pair, other, effect, and weight, and calculates a genetic risk score applying the formula $\sum \beta_{i}G_{i}$, where $\beta_i$ is the weight for snp $i$ and $G_i$ is the individuals's imputed genotype for SNP $i$. This can be run using

`a=generate_grs('Example_GRS')`.

Note that this is the same file used in the extract_snp_bulk function. This is intentional, as it allows users to check each any individual SNP in the GRS easily without generating a new file.

## For internal Exeter use
