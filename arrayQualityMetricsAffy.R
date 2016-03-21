setwd("~/a/data")
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("arrayQualityMetrics")

library(nugohs1a520180cdf)
library(nugohs1a520180.db)


install.packages(pkgs = c("nugohs1a520180hsentrezg.db", "nugohs1a520180hsentrezgcdf", "nugohs1a520180hsentrezgprobe"), repos = "http://nmg-r.bioinformatics.nl/bioc")

library(nugohs1a520180hsentrezgcdf)
library(nugohs1a520180hsentrezg.db)
#select(hgu133plus2.db, yourprobeids, c("SYMBOL","ENTREZID"))                 

gseString = "GSE53291"
gseUntaredFolder = paste(gseString,  "_untared", sep = "")
celFilesPath = paste("~/a/data", gseUntaredFolder, sep = "/")
setwd(celFilesPath)

library(affy)
library(sva)

pd = read.AnnotatedDataFrame(filename="pheno.txt")
affyBatch = affy::ReadAffy(phenoData=pd,
                          sampleNames=pd$SampleAccessionNumber)
                        #  cdfname = "nugohs1a520180_hs_entrezg")

cel_table = read.table("GSM12890
                       02_Placenta_CON_female_1.CEL")

# according to ncbi the platform has 23941 probesets
# according to ncbi the platform for custom cdf has 18525 probesets
# following values are "standard/custom cdf"
length(featureNames(affyBatch))                 # 23941/17451 probesets 
nrow(affyBatch@assayData$exprs)                 # 535824/535824 probes

affyBatch.rma <- affy::rma(affyBatch)

affyBatch.rma@featureData@data                  # 23941/17451 probesets 
length(rownames(affyBatch.rma@assayData$exprs)) # 23941/17451 probesets 
nrow(affyBatch.rma@assayData$exprs)             # 23941/17451 probesets 


library(AnnotationDbi)
probesetsID<-rownames(affyBatch.rma@assayData$exprs)

# custom cdf
length(probesetsID)
length(keys(nugohs1a520180hsentrezg.db, "SYMBOL")) #  59476
length(keys(nugohs1a520180hsentrezg.db, "ENTREZID")) # 17349
probesetsID_EntrezID<-select(nugohs1a520180hsentrezg.db, probesetsID, "ENTREZID")
length(probesetsID_EntrezID$PROBEID) # amount of annotated probesets in mapping 17451

# standard cdf
length(probesetsID)
length(keys(nugohs1a520180.db, "SYMBOL")) # 59476
length(keys(nugohs1a520180.db, "ENTREZID")) # 16922
probesetsID_EntrezID<-select(nugohs1a520180.db, probesetsID, "ENTREZID")
length(probesetsID_EntrezID$PROBEID) # amount of annotated probesets in mapping 24026


design <- model.matrix(~0 + colnames(exprs))
library(limma)
fit <- lmFit(exprs(exprs), design)
ebFit <- eBayes(fit)

# exclude batch effect using day of experiment as a batch
# pd$ScanDate <- substr(exprs@protocolData$ScanDate, 1, 10)
# batch = pd$ScanDate
# mod = model.matrix(~as.factor(SampleType), data=pData(pd))
# combat_edata = sva::ComBat(dat=exprs(exprs), 
#                       batch=batch,
#                       mod=mod,
#                       par.prior=TRUE,
#                       prior.plots=FALSE)
# exprs(exprs) <- combat_edata

arrayQualityMetrics::arrayQualityMetrics(expressionset = affyBatch.rma,
                                         outdir = "QC_affy_custom_cdf",
                                         force = TRUE,
                                         do.logtransform = FALSE,
                                         intgroup = "Target")

