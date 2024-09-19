## Makes mathatten plot from filtered vcf file and meta data
#compares chromosome region and fst

library(vcfR)
library(tidyverse)
library(qqman)

X11.options(type="cairo")

vcf <- read.vcfR("~/Rprojects/eco_genomics/eco_genomics/outputs/vcf_final.filtered.vcf.gz")

### metadata
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

## check if the name number of rows... (oh no!! they don't!!)
vcf
dim(meta)

#let's fix it! (make meta have only individuals that match!)
meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]
head(vcf)


## calculate diversity statistics using genetic_diff
vcf.div <- genetic_diff(vcf, 
                        pops=as.factor(meta2$region),
                        method="nei")

#take first 8 unique chromosomes (no scaffolds or repeats!)
chr.main <- unique(vcf.div$CHROM)[1:8]
#assign a number 
chrnum <- as.data.frame(cbind(chr.main, seq(1,8, 1)))

## reorganize data for manhatten blot (join)
vcf.div.MHP <- left_join(chrnum, vcf.div, join_by(chr.main==CHROM))

## add a blank SNP line
vcf.div.MHP <- vcf.div.MHP %>% filter(Gst>0) %>% mutate(SNP=paste0(chr.main,"_",POS))

#make strings numeric
vcf.div.MHP$V2 = as.numeric(vcf.div.MHP$V2)
vcf.div.MHP$POS = as.numeric(vcf.div.MHP$POS)

### MAKE manhattan plot!!!! yipeee
manhattan(vcf.div.MHP,
          chr="V2",
          bp="POS",
          p="Gst",
          col=c("lightblue4","orange3"),
          logp=F, 
          ylab="Fst among regions",
          suggestiveline=quantile(vcf.div.MHP$Gst, 0.999))
