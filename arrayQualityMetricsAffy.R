setwd("~/a/data")
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("arrayQualityMetrics")
biocLite("nugohs1a520180hsentrezg.db")
library(nugohs1a520180cdf)
library(nugohs1a520180hsentrezg.db)

biocLite("sva")

gseString = "GSE53291"


gseUntaredFolder = paste(gseString,  "_untared", sep = "")
celFilesPath = paste("~/a/data", gseUntaredFolder, sep = "/")
setwd(celFilesPath)

library(affy)
library(sva)
pd = read.AnnotatedDataFrame(filename="pheno.txt")
affyData = affy::ReadAffy(phenoData=pd,
                          sampleNames=pd$SampleAccessionNumber,
                          cdfname = "NuGO_Hs1a520180")
head(featureNames(affyData)) # probe_ids


nexprs <- affy::rma(affyData)

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

arrayQualityMetrics::arrayQualityMetrics(expressionset = exprs,
                                         outdir = "QC_affy_custom_cdf",
                                         force = TRUE,
                                         do.logtransform = FALSE)
