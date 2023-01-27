#------------------------------------------------------------------------------
# Script Name: extract_snp_bulk
# Purpose: A function to extract multiple SNPs from DNA Nexus' imputed genotype data using a file input
# Author: Dr. Harry Green, University of Exeter
# Date Created: 17/01/23
# Dependencies: depends on the package rbgen being installed. Use install.packages( "http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.5.tgz", repos = NULL, type = "source" ) then library('rbgen') to ensure the correct version is installed
# Notes: file_in should be a tab separated file with the columns chromosome, bp, other, effect, weight, where other and effect are the allele codes. If the names are not consistent, they will be renamed, it's the order that matters
#------------------------------------------------------------------------------

extract_snp_bulk=function(file_in){
  
  
  # This function just creates a list of file names of bgen files in UKBB. This should be constant across all projects, I think
  create_filenames=function(){
    return(c(
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c1_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c2_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c3_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c4_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c5_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c6_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c7_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c8_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c9_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c10_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c11_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c12_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c13_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c14_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c15_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c16_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c17_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c18_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c19_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c20_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c21_b0_v3.bgen",
      "../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c22_b0_v3.bgen"))
  }
  filenames=create_filenames()
  
  sample=read.table("../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c1_b0_v3.sample",header=T) #this file has all the sample ids in it, but there's a dummy line at the start
  eid=sample$ID_1[2:nrow(sample)]
  
  snps_in=read.table(file_in,header=T)
  snps_in=snps_in[,c(1,2)]
  
  #people should be putting these headers on anyway, but in case they're wrong, I've relabelled them
  names(snps_in)=c('chr','bp')
  nsnps=nrow(snps_in)
  
  #
  genotypes=data.frame(eid=eid)
  variants=data.frame(chromosome=double(),position=double(),rsid=character(),number_of_alleles=character(),allele0=character(),allele1=character())
  
  for (i in 1:22){
    print(paste('extracting SNPs on chromosome',i))
    
    snps_in_i=snps_in[snps_in$chr==i,] #trying not to use dplyr, but this is just the snps on chromosome i
    
    if (nrow(snps_in_i)>0){ #if there's nothing on the chromosome just skip it
      #the rbgen package wants the chromosome in two digit form, e.g. '08'
      if (i<10){
        chr=paste('0',i,sep='')
      }else{
        chr=as.character(i)
      }
      
      #this tells rbgen where to look
      ranges = data.frame(
        chromosome = rep(chr,nrow(snps_in_i)),
        start = snps_in_i$bp,
        end = snps_in_i$bp
      )
      data=bgen.load(filenames[i], ranges )# this pulls out the data for all snps on the chromosome. It has to be by chromosome because the dna nexus data is stored in one file per chromosome
	  
	  variants=rbind(as.data.frame(variants),as.data.frame(data$variants))
      
      for (j in 1:nrow(snps_in_i)){
		colname=paste(i,':',snps_in_i$bp[j],sep='')
        #if it doesn't find a variant it causes problems, so I need to match the base pair
        datavar=which(snps_in_i$bp[j]==data$variants$position) #this is the row in the extracted data that corresponds to the variant j in grs_chr_i
        
        if (length(datavar>0)){ #so only if there's a matching base pair
			for (k in datavar){
				 paste(i,':',data$variants$position[k],'_',data$variants$allele0[k],'_',data$variants$allele1[k],sep='')
				 mat=data$data[k,,]
				 geno=as.numeric(mat[,2]+2*mat[,3])
				 genotypes[[colname]]=geno
				 }
        }
      }
    }
  }
  
  out=list(genotypes=genotypes,variants=variants)
  return(out)
}




