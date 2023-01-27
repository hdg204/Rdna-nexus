install.packages( "http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.5.tgz", repos = NULL, type = "source" )
library('rbgen')
library('dplyr')

source_url("https://raw.githubusercontent.com/hdg204/Rdna-nexus/main/extract_snp.R") 
source_url("https://raw.githubusercontent.com/hdg204/Rdna-nexus/main/extract_snp_bulk.R") 
source_url("https://raw.githubusercontent.com/hdg204/Rdna-nexus/main/generate_grs.R") 
system("curl -o Example_GRS https://raw.githubusercontent.com/hdg204/Rdna-nexus/main/Example_GRS")
