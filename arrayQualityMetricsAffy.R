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
                          sampleNames=pd$SampleAccessionNumber,
                          cdfname = "nugohs1a520180_hs_entrezg")

cel_table = read.table("GSM1289002_Placenta_CON_female_1.CEL")
<<<<<<< HEAD
nrow(cel_table)   # 21933 rows in cel file=probes in cel file?
=======
nrow(cel_table)   # 21933 rows in cel file?
>>>>>>> cef4df7ab239bb8f92e21fc65b72d0c37e2ae446

# according to ncbi the platform has 18525 probesets
length(featureNames(affyBatch))                 # 17451 named probesets amount?
nrow(affyBatch@assayData$exprs)                 # 535824 probes amount?

affyBatch.rma <- affy::rma(affyBatch)

affyBatch.rma@featureData@data                  # probesets amnt 17451
length(rownames(affyBatch.rma@assayData$exprs)) # unique probeset names amnt 17451
nrow(affyBatch.rma@assayData$exprs)             # probeset amnt 17451


library(AnnotationDbi)
probesetsID<-rownames(exprs@assayData$exprs)
<<<<<<< HEAD
length(keys(nugohs1a520180hsentrezg.db, "SYMBOL")) # amount of probesets in annotation library or 59476
=======
length(probesetsID)
length(contents(nugohs1a520180hsentrezgSYMBOL)) # amount of probesets in annotation library 17451
length(keys(nugohs1a520180hsentrezg.db, "SYMBOL")) # or 59476
>>>>>>> cef4df7ab239bb8f92e21fc65b72d0c37e2ae446
length(keys(nugohs1a520180hsentrezg.db, "ENTREZID")) # or 17349, depends on mapping
probesetsID_EntrezID<-select(nugohs1a520180hsentrezg.db, probesetsID, "ENTREZID")
length(probesetsID_EntrezID$PROBEID) # amount of annotated probesets in microarray 17451



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

