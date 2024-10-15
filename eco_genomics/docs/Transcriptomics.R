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
percentVar <- round(100*attr(pcaData, "percentVar"))

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





########################
## DAY 2 TRANSCRIPTOMICS
########################


library(pheatmap)
resultsNames(dds)
##[1] "Intercept" "DevTemp_D22_vs_D18" "FinalTemp_A33_vs_A28" "FinalTemp_BASE_vs_A28"

#   :P

res_D22vsD18 <- results(dds, name="DevTemp_D22_vs_D18", alpha=0.05)

#orders by most signifcant differential expression
res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$pvalue),]

head(res_D22vsD18)


#forkle
#sporkle
#crumb
#leefie
#fhawks

d <- plotCounts(dds, gene="TRINITY_DN147593_c1_g4_i1", int=(c("DevTemp", "FinalTemp")), returnData=T)
d

p <- ggplot(d, aes(x=DevTemp, y=count, color=DevTemp, shape=FinalTemp)) + 
  theme_minimal() + theme(text=element_text(size=10),panel.grid.major=element_line(colour = "slateblue4")) + 
  geom_point(position = position_jitter(w=0.2,h=0),size=3)
p

p.MA <- plotMA(res_D22vsD18, ylim=c(-4,4))

#fourtnite!!!!!!!!!!!!!!!!!!!!

#gadagada ee gadagada oh
   ____
 _(
(__)

#volcano plot
res_df <- as.data.frame(res_D22vsD18)

res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")
res_df

p2 <- ggplot(res_df, aes(x= log2FoldChange, y=-log10(padj), color = Significant))+
  geom_point(alpha=0.8 #,position = position_jitter(w=0.2,h=.2) #goofy ahh jitter
             )+
  scale_color_manual(values = c("slateblue","hotpink")) + 
  labs(x="Log2 Fold Change", y="log10 Adjusted P Value", title = "Volcano Plot") +
  theme(legend.position = "top") +
  geom_hline(yintercept=-log10(0.05), linetype = "dashed", color = "green3") +
  geom_vline(xintercept = c(-1,1), linetype="dashed", color = "green3") #+ ylim(0.2,10)
p2


vsd <- vst(dds, blind=F)
topgenes <- head(rownames(res_D22vsD18), 20)
mat <- assay(vsd)[topgenes,]
df <- as.data.frame(colData(dds)[,c("DevTemp", "FinalTemp")])
pheatmap(mat, annotation_col=df, show_rownames=F, cluster_cols=T, cluster_rows=T)
