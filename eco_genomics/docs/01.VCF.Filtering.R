library(vcfR)
#hi vcfr!! why are you working??? I love you fr fr <3
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

getwd() #basically pwd
list.files() #basically ls

#read vcf files
vcf <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/Centaurea_filtered.vcf.gz")
head(vcf)

#get reference genome
dna <- ape::read.dna("/gpfs1/cl/pbio3990/PopulationGenomics/reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")

#annotation
gff <- read.table("/gpfs1/cl/pbio3990/PopulationGenomics/reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")


#Put them all together for chromosome 1
chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna, ann=gff)
###Makes a plot!!1
plot(chr1)

chromoqc(chr1, xlim=c(1e1, 1.1e8))
