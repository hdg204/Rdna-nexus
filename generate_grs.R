#------------------------------------------------------------------------------
# Script Name: generate_grs
# Purpose: A function to calculate a genetic risk score on the DNA Nexus platform from the RStudio Workbench implementation
# Author: Dr. Harry Green, University of Exeter
# Date Created: 17/01/23
# Dependencies: depends on the package rbgen being installed. Use install.packages( "http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.5.tgz", repos = NULL, type = "source" ) then library('rbgen') to ensure the correct version is installed
# Notes: file_in should be a tab separated file with the columns chromosome, bp, other, effect, weight, where other and effect are the allele codes. If the names are not consistent, they will be renamed, it's the order that matters
#------------------------------------------------------------------------------

generate_grs=function(file_in){
  
  dosage=extract_snp_bulk(file_in)
  dosage_matrix=dosage$genotypes[,2:ncol(dosage$genotypes)]
  
  sample=read.table("../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c1_b0_v3.sample",header=T) #this file has all the sample ids in it, but there's a dummy line at the start
  eid=sample$ID_1[2:nrow(sample)]
  
  grs_snps=read.table(file_in,header=T)
  grs_snps=grs_snps%>%arrange(chr,bp)
  
  #some of the snps will have not been read in by the UK Biobank, so I need to make sure that only the SNPs in dosage_matrix are in grs_snps. The thing is the dosage matrix has snps ordered by the bgen file, and the snp input can be either way around
  dosage_snp_names=names(dosage_matrix)
  input_snp_names=paste(grs_snps$chr,':',grs_snps$bp,'_',grs_snps$other,'_',grs_snps$effect,sep='')
  
  #some will match, some will match if a flip is made
  matches=rep(NA,length(input_snp_names)) #for each snp in the input file, index is its position in the output file
  flip_matches=rep(NA,length(dosage_snp_names)) #for each snp in the dosage, flip is if a flip was required to match it to the input file
  
  for (i in 1:length(input_snp_names)){
	for(j in 1:length(dosage_snp_names)){
		if (input_snp_names[i]==dosage_snp_names[j]){ #first try a straight match
			matches[i]=j
		}else{ #if they don't match, I'll try a flip.
			flipname=paste(unlist(strsplit(input_snp_names[i],'_')[[1]][c(1,3,2)]),collapse='_')
			if (flipname==dosage_snp_names[j]){
				matches[i]=j
				flip_matches[j]=j
			}
		}	
	}
  }
  
  missing_snps=input_snp_names[is.na(matches)]
  
  #trim the input file to only the ones that are there. We don't need the ones that didn't match now I have them in their own output variable
  grs_snps=grs_snps[!is.na(matches),]
  
  #for flipped snps, I will flip the dosage file by 2- instead of flipping the weights, because then it will work even if odds ratios are given instead of log-odds BUT IT DOESN'T WORK FOR SEX CHROMOSOMES
  # dosage_matrix[,!is.na(flip_matches)]=2-dosage_matrix[,!is.na(flip_matches)]
  
  grs_snps$weight[!is.na(flip_matches)]=-grs_snps$weight[!is.na(flip_matches)]
  
  a=as.matrix(dosage_matrix)
  b=matrix(grs_snps$weight) 
  grs=a%*%b #The entire GRS is made by this neat matrix multiplication, where the dosage table is multiplied by the vector of weights to give a vector of risk scores
  grs_df=data.frame(eid=eid,grs=grs)
  list(grs=grs_df,missing=missing_snps,snp_data=dosage)
  return(grs)
}
