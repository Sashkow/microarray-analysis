setwd("~/a/data/GSE30186_untared/")
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("lumi")
biocLite("illuminaHumanv4.db")

library(lumi)
library(illuminaHumanv4.db)
fileName = 'GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt'
fileName1 = "GSE30186_non_normalized.txt"

lumiBatch = lumiR.batch(fileName1, sampleInfoFile = "sampleInfo.txt")

length(lumiBatch@assayData$exprs)

limiProcessedBatch = lumiExpresso(lumiBatch)



arrayQualityMetrics::arrayQualityMetrics(expressionset = limiProcessedBatch,
                                         outdir = "QC_lumi_processesd",
                                         force = TRUE,
                                         do.logtransform = FALSE,
                                         intgroup = "Target")


summary(limiProcessedBatch, "QC")


lumi_probesetid_list = featureNames(limiProcessedBatch)

probeid_to_gene_mapping = select(illuminaHumanv4.db, lumi_probesetid_list, "ENTREZID")



# drawing barplot: amount of probesets per gene -- amount of genes that have such frequency
library(ggplot2)
View(probeid_to_gene_mapping)
gene_frequency_table = table(probeid_to_gene_mapping$ENTREZID)
# genes of all probeset amounts 
gene_frequency_table = as.data.frame(gene_frequency_table)
# genes of probeset amounts > 10
filtered_gene_frequency_table = data.frame()
for (i in 1:nrow(gene_frequency_table)){
  row <- gene_frequency_table[i,]
  if (row$Freq>10){
    filtered_gene_frequency_table = rbind(filtered_gene_frequency_table, row)
  }
}

c <- ggplot(filtered_gene_frequency_table, aes(factor(Freq)))
c <- ggplot(gene_frequency_table, aes(factor(Freq)))
c + geom_bar()

lumiBatchFrame <- data.frame(lumiBatchN)
annot <- data.frame(SYMBOL=sapply(contents(illuminaHumanv4SYMBOL), paste, collapse=", "))


all <- merge(annot, lumiBatchFrame, by.x=0, by.y=0, all=T)


design <- model.matrix(~0 + colnames(lumiBatchN))
View(design)
library(limma)
fit <- lmFit(lumiBatchN, design)
ebFit <- eBayes(fit)








