extract_snp=function(chr,bp){
  # This function just creates a list of file names of bgen files in UKBB. This should be constant across all projects, I think
  bgen_file=paste("../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c",chr,"_b0_v3.bgen",sep='')
 
  sample=read.table("../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_c1_b0_v3.sample",header=T) #this file has all the samples in it, but there's a dummy line at the start
  eid=sample$ID_1[2:nrow(sample)]
  
  #the rbgen package wants the chromosome in two digit form, e.g. '08'
  if (chr<10){
    chr_str=paste('0',chr,sep='')
  }else{
    chr_str=as.character(chr)
  }
  
  # SEX CHROMOSOMES BREAK THINGS. I fix.
  if (chr==23){
	# Some people aren't in the sex chromosome bgen file so using the sample file from chromosome 1 doesn't work
	sample=read.table("../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_cX_b0_v3.sample",header=T)
	bgen_file="../../mnt/project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22828_cX_b0_v3.bgen"
	chr_str='X'
	eid=sample$ID_1[2:nrow(sample)]
  }
  
  #this tells rbgen where to look
  ranges = data.frame(
    chromosome = chr_str,
    start = bp,
    end = bp
  )
  data=bgen.load(bgen_file, ranges )
  
  if (nrow(data$variants)==0){
	return('SNP not found')
  }else{
  mat=data$data[1,,]
  genotypes=as.numeric(mat[,2]+2*mat[,3])
  
  out=list(snp_data=data.frame(eid,genotypes),variant_info=data$variants)
  return(out)
  }
}
