# Rdna-nexus

## Installation

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

## extract_snp_bulk

While many SNPs can be extracted using `extract_snp`, this is slow as it would need to read the bgen for each SNP. The most efficient way is to read the bgen for each chromosome, since these are split in the UKBB data. This function is built for extracting multiple SNPs from a file.

The file should be a tab separated file, with the first two columns being chromosome and base pair position. An example of this is given in https://github.com/hdg204/Rdna-nexus/blob/main/Example_GRS, where only the first two columns are used.

After installation, this function can be run using:

`a=extract_snp_bulk('Example_GRS')`

It can be a little slow but it prints to screen which chromosome it's on. The output is similar to extract_snp, except there are multiple columns for the genotypes, named 'chr:bp_A0_A1', and multiple rows for the SNP info, one for each SNP.

## generate_grs

## For internal Exeter use
