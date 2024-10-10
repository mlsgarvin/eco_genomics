## analyze RNA seq using DESeq2

library(DESeq2)
library(tidyverse)

setwd("~/Rprojects/eco_genomics/eco_genomics/docs")

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt", header=T, row.names=1)

tail(countsTable)

##make everything a whole number
countsTableRound <- round(countsTable)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt", header = T, stringsAsFactors = T, row.names = 1)


#make some plots
barplot(colSums(countsTableRound), names.args=colnames(countsTableRound), 
        cex.names=0.5, las=2, ylim=c(0,30000000))

abline(h=mean(colSums(countsTableRound)), col="lightgreen", lwd=2)

##lock in right now yalllll
## feelin vacant... yeaaaaa
## when i wanted your love i wasnt strong enough
### exolore counts matrix!!

colSums(countsTableRound)
mean(colSums(countsTableRound))
## aim for 20 million at design, and after filtering we have 18 million!! great work!!

barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound),
        cex.names=0.5, las = 2, ylim=c(0,30000000))



## explore da matrix a little bit moar

## average number of counts per gene
mean(rowSums(countsTableRound))
median(rowSums(countsTableRound))

#mean across columns
apply(countsTableRound,2,mean)




#on the pavement
#talk to me
#heart's feelin vacant

#coding in poems
#poems in code
# where would I be
#without seeing you both

#my love us eternal
#but my feelings are not
#give what you can
#cuz it's all that youve got

#marshmellow skies
#and a bright purple lake
# going for a swim
#after our last slice of cake



dds <- DESeqDataSetFromMatrix(countData= countsTableRound, 
                              colData=conds, design = ~ DevTemp + FinalTemp)

dds <- dds[rowSums(counts(dds) >= 10) >= 15,]
nrow(dds)


#run da model!!!
dds <- DESeq(dds)

###List results?!?!?!
resultsNames(dds)

#[1] "Intercept"             "DevTemp_D22_vs_D18"    "FinalTemp_A33_vs_A28" 
#[4] "FinalTemp_BASE_vs_A28"

##PCA plot time (yipeee!!!!)
# step 1: transform
vsd <- #variance stabilized data
  vst(dds,blind=F)

pcaData <- plotPCA(vsd, intgroup=c("DevTemp", "FinalTemp"), returnData=T)
percentVar <- round(100*attr(pcsData, "percentVar"))

final_temp_colors <- c("BASE" = "grey15", "A28" = "hotpink", "A33" = "red3")
shapes_choose <- c("D18" = 16, "D22" = 18)

p <- ggplot(pcaData, aes(PC1, PC2, color=FinalTemp, shape=DevTemp)) + 
  geom_point(size = 5) +
  scale_shape_manual(values = shapes_choose) +
  scale_color_manual(values = final_temp_colors) +
  labs(x =paste0('PC1: ', percentVar[1], ' %'), y = paste0('PC2: ', percentVar[2], '% ')) +
  theme_bw(base_size = 16)

p

## Fish in the pool
#by Miles Garvin 

#a fish
#in the pool

#how did he get there?
#who let him in?

#what are you doing?
#little fish
#in the pool


#the vaccum rumbles
#and bubbles rise from the water's edge
#and you swim
#swim 
#swim

#not for much longer

#how will you get out?
#why wont you let me help?
#come here, little fishy
#don't you want out?

#here lies a fish
#at the bottom of a pool
#how did he get there?
#in a world ever cruel