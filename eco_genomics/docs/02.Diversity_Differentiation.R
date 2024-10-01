## Makes mathatten plot from filtered vcf file and meta data
#compares chromosome region and fst
#9/19

library(vcfR)
library(tidyverse)
library(qqman)

X11.options(type="cairo")
options(bitmapType="cairo")

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


write.vcf( )





### Run admixture analysis and create plots
library(LEA)
CentAdmix <- snmf("vcf_final.filtered.vcf.gz", #### change to be the right one
                  K=1:10,
                  entropy=T,
                  repetitions= 5,
                  project="new") ##if ur adding to this later, change project to continue

par(mfrow=c(2,1))
plot(CentAdmix)
plot(CentAdmix$eigenvalues[1:10],
     ylab="Eigen Values",
     xlab="Number of PCs",
     col="lightblue4",
     main="PCA")

dev.off() ### undoes parameters (yippeeeeeeeeee)


myK=4

CE = cross.entropy(CentAdmix, K=myK)
best = which.min(CE) #chooses the LOWEST run --> best



myKQ = Q(CentAdmix, K=myK, run=best)
myKQOmega = cbind(myKQ, meta2)

my.colors = c("cyan", "tomato", "forestgreen", "violet", "goldenrod","purple3", "grey23", "blue2", "lightgreen", "maroon4", "darkred", "orange3", "yellow1", "white", "grey4", "tan", "tan4")
myKQmeta = as_tibble(myKQmeta) %>% 
  group_by(continent) %>% 
  arrange(region, pop, .by_group = FALSE)

barplot(as.matrix(t(myKQmeta[,1:myK])),
        border=NA,
        space=0,
        col=my.colors[1:myK],
        xlab="Geographic Reigons",
        ylab="Ancestry Proportions",
        main=paste0("Ancestry Matrix K=",myK))
axis(1,
     at=1:length(myKQmeta$region),
     labels=myKQmeta$region,
     tick=F,
     cex.axis=0.5,
     las=3)

barplot(my.colors,
        border=NA,
        space=0.5,
        col=my.colors[1:length(my.colors)],
        xlab="Geographic Reigons",
        ylab="Ancestry Proportions",
        main=paste0("Ancestry Matrix K=",myK))




