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

lumiProcessedBatch = lumiExpresso(lumiBatch)

arrayQualityMetrics::arrayQualityMetrics(expressionset = limiProcessedBatch,
                                         outdir = "QC_lumi_processesd",
                                         force = TRUE,
                                         do.logtransform = FALSE,
                                         intgroup = "Target")
summary(limiProcessedBatch, "QC")

lumi_probesetid_list = featureNames(lumiProcessedBatch)
length(lumi_probesetid_list) # 47318

probeid_to_gene_mapping = select(illuminaHumanv4.db, lumi_probesetid_list, "ENTREZID") 
length(probeid_to_gene_mapping$PROBEID) # 50627 != 47318
length(unique(probeid_to_gene_mapping$PROBEID)) # 47318



# exclude probesets that do not map on genes
probeid_to_gene_mapping = probeid_to_gene_mapping[!is.na(probeid_to_gene_mapping$ENTREZID),] # 39191
# exclude probestts that map on multiple gnes
duplicated_probeid = unique(probeid_to_gene_mapping[duplicated(probeid_to_gene_mapping$PROBEID),]$PROBEID) # amount of unique probeids that have duplicaties 1726
probeid_to_gene_mapping = probeid_to_gene_mapping[!(probeid_to_gene_mapping$PROBEID %in% duplicated_probeid),] # 34156
# separate probesets that map on unique gene and those that map on multiple genes
unique_entrezid = unique(probeid_to_gene_mapping[duplicated(probeid_to_gene_mapping$ENTREZID),]$ENTREZID) # 7703
probeid_to_gene_mapping_duponly = probeid_to_gene_mapping[probeid_to_gene_mapping$ENTREZID %in% unique_entrezid,] # 20818
probeid_to_gene_mapping_no_dupl = probeid_to_gene_mapping[!(probeid_to_gene_mapping$ENTREZID %in% unique_entrezid),] # 13338


lumiBatchFrame = data.frame(lumiProcessedBatch@assayData$exprs)
geneToExpressionValueMappingDuplicate = lumiBatchFrame[0,] # empty data frame with same column
geneToExpressionValueMappingNoDuplicate = data.frame()
lumiBatchFrame$PROBEID <- rownames(lumiBatchFrame)

# merge probeid_to_gene_mapping_no_dupl and lumiBatchFrame by probesetid
ha = probeid_to_gene_mapping_no_dupl
hb = lumiBatchFrame
geneToExpressionValueMappingNoDuplicate = merge(ha, hb)
geneToExpressionValueMappingNoDuplicate$PROBEID = NULL
length(colnames(geneToExpressionValueMappingNoDuplicate))


for (gene_name in unique_entrezid){
  probesetidList = probeid_to_gene_mapping_duponly[probeid_to_gene_mapping_duponly$ENTREZID == gene_name,]$PROBEID
  expressonValuesForGene = lumiBatchFrame[lumiBatchFrame$PROBEID %in% probesetidList,]
  apply(expressonValuesForGene[1:ncol(expressonValuesForGene)-1],2, as.numeric)
  # labeled vector of average expression values for probesets pointing go a gene for each data set
  a = apply(expressonValuesForGene[1:ncol(expressonValuesForGene)-1],2, mean)
  geneToExpressionValueMappingDuplicate[nrow(geneToExpressionValueMappingDuplicate)+1,] <- a
}

geneToExpressionValueMappingDuplicate$ENTREZID = unique_entrezid

length(colnames(geneToExpressionValueMappingDuplicate))

geneToExpressionValueMapping = rbind(geneToExpressionValueMappingNoDuplicate,geneToExpressionValueMappingDuplicate)

nrow(geneToExpressionValueMapping) # 21041
# ANSWER
View(geneToExpressionValueMapping) 

#test ENTREZID is unique in geneToExpressionValueMapping
duplicate_entrezid = unique(geneToExpressionValueMapping[duplicated(geneToExpressionValueMapping$ENTREZID),]$ENTREZID)
length(duplicate_entrezid) == 0



# drawing barplot: amount of probesets per gene -- amount of genes that have such frequency
library(ggplot2)
View(probeid_to_gene_mapping)
gene_frequency_table = table(probeid_to_gene_mapping$ENTREZID)
# genes of all probeset amounts 
gene_frequency_table = as.data.frame(gene_frequency_table)
# genes of probeset amounts > 10 (to make sence of scale)
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

design <- model.matrix(~0 + colnames(lumiBatchN))
View(design)
library(limma)
fit <- lmFit(lumiBatchN, design)
ebFit <- eBayes(fit)








