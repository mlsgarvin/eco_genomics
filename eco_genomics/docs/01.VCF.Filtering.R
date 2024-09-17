library(vcfR)
library(tidyverse)
#hi vcfr!! why are you working??? I love you fr fr <3
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

getwd() #basically pwd
list.files() #basically ls

#read vcf files
vcf <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/Centaurea_filtered.vcf.gz")
vcf
head(vcf, 10)

#get reference genome
dna <- ape::read.dna("/gpfs1/cl/pbio3990/PopulationGenomics/reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")

#annotation
gff <- read.table("/gpfs1/cl/pbio3990/PopulationGenomics/reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")


#Put them all together for chromosome 1
chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna, ann=gff)
###Makes a plot!!1
plot(chr1)

#open for writing to a pdf
pdf(file="~/Rprojects/eco_genomics/eco_genomics/figures/Chromoplot.pdf")

chromoqc(chr1, xlim=c(1e1, 1.1e8)) #then run command

#closes pdf writer
dev.off()




##### Day 2
####
####

DP <- extract.gt(vcf, element="DP", as.numeric=T)
dim(DP)

#show first 5 rows, 10 individuals
DP[1:5,1:10]

quantile(DP)

##assign NA to values that are 0 (missing)
DP[DP==0] <- NA

quantile(DP, na.rm=T)


### Visualize missingness in VCF file
heatmap.bp(DP)   ### Not working :(

##subset. use only first 1000 loci
heatmap.bp(DP[1:1000,], rlabels=F, clabels=F)


##filter by DEPTH
library(SNPfiltR)
#low depth, filters out genotypes with a depth < 3
vcf.filt <- hard_filter(vcf, depth=3)

#high depth filters out genotypes with a depth > 60
vcf.filt <- max_depth(vcf.filt, maxdepth=60)


##Fileter by MISSINGNESS
meta <- read.csv("metadata/meta4vcf.csv", header=T)
head(meta2, 7)

## take all rows, only colums 1 and 4 (Concatonated)
meta2 <- meta[,c(1,4)]
#relabel columns
names(meta2) <- c("id","pop")

#make them into factors?? troubleshooting
meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)


## filter out INDIVIDUALS with > 75% missing data
vcf.filt.indMiss <- missing_by_sample(vcf.filt,
                                      popmap=meta2,
                                      cutoff=0.75)
#filter to be bialellic
vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)

#get rid of NAs
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)

##snp-wise missingness filter
vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff=0.5)


## save final file!!
write.vcf(vcf.filt.indSNPMiss,
          "~/Rprojects/eco_genomics/eco_genomics/outputs/vcf_final.filtered.vcf.gz")
