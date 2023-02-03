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


sample=read.table("../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c1_b0_v3.sample",header=T) #this file has all the sample ids in it, but there's a dummy line at the start
eid=sample$ID_1[2:nrow(sample)]

grs_snps=read.table(file_in,header=T)
grs_snps=grs_snps%>%arrange(chr,bp)

nsnps=nrow(grs_snps)

missing=rep(0,nsnps)
flip=rep(0,nsnps)
multi=rep(0,nsnps)

#pulling out chromosome from dosage is a bit difficult due to stringiness.
chrs=as.numeric(dosage$variants$chromosome)
chrs[dosage$variants$chromosome=='X']=23

for (i in 1:nsnps){
#for each snp, m1 will be 1 if the variants and alleles match, and 0 if it doesn't
#for each snp, m2 will be 1 if the variants match alleles but alleles are swapped, and 0 if it doesn't
# m1 and m2 will both be 0 if there is no match even after a flip
# m3 is how many chr and bp matches there are. >1 for multi-allelic SNPs
	m1=(sum(
			  grs_snps$chr[i]==chrs&
			  grs_snps$bp[i]==dosage$variants$position&
			  grs_snps$effect[i]==dosage$variants$allele1&
			  grs_snps$other[i]==dosage$variants$allele0
			  )
	)
	m2=(sum(
			  grs_snps$chr[i]==chrs&
			  grs_snps$bp[i]==dosage$variants$position&
			  grs_snps$effect[i]==dosage$variants$allele0&
			  grs_snps$other[i]==dosage$variants$allele1
			  )
	)
	m3=(sum(
			  grs_snps$chr[i]==chrs&
			  grs_snps$bp[i]==dosage$variants$position
			  )
	)
	if (m1+m2==0){
		missing[i]=1
	}
	if (m2==1){
		flip[i]=1
	}
	if (m3>1){
		multi[i]=1
	}
}

missing_df=grs_snps[which(missing==T),]

# missing, flip, and multi refer to rows within grs_snps. If grs_snps changes size, these variables will no longer be accurate.

#I will first flip anything that needs flipping, because this will not change the number of rows in the table
grs_snps$weight[flip==1]=-grs_snps$weight[flip==1]

# if there are multi allelic snps, these need to be removed from the DOSAGE file but NOT from the snp list. There will be 2 columns for that SNP in the dosage file and I need to delete the correct one

# this code looks similar to the other one, but the referencing is going the other way

rem=c()
for (i in which(multi==1)){
	# just to shorten code, store the correct chr bp oth eff into variables
	ch=grs_snps$chr[i]
	bp=grs_snps$bp[i]
	oth=grs_snps$other[i]
	eff=grs_snps$effect[i]
	chrbpmatch=chrs==ch & dosage$variants$position==bp #is it the right coordinate?
	allelematch=dosage$variants$allele1==eff & dosage$variants$allele0==oth #are the alleles correct?
	flipmatch=dosage$variants$allele0==eff & dosage$variants$allele1==oth #are the alleles flipped?
	
	wrongallele=which(chrbpmatch & !(allelematch | flipmatch)) #if the chr matches, but the alleles don't even after a flip, we're looking at the wrong allele
	# remove it from the variants list
	rem=c(rem,wrongallele)
}
dosage$variants=dosage$variants[-rem,]
dosage$genotypes=dosage$genotypes[,-(rem+1)]

# now I will remove the missing ones. The length of the dosage file should be sum(missing) less than the length of the snp file
grs_snps=grs_snps[-which(missing==T),]

a=as.matrix(dosage$genotypes[,2:ncol(dosage$genotypes)])
b=as.matrix(grs_snps$weight) 
grs=a%*%b

grs=data.frame(eid=eid,grs=grs)%>%head()

return(out=list(grs=grs,dosage=dosage,missing_snps=missing_df))
}


