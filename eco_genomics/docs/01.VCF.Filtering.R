library(vcfR)
#hi vcfr!! why are you working??? I love you fr fr <3
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

getwd() #basically pwd
list.files() #basically ls

#read vcf files
vcf <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/Centaurea_filtered.vcf.gz")
head(vcf)
