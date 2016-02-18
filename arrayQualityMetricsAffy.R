
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("arrayQualityMetrics")
biocLite("sva")
gseString = "GSE53291"

gseUntaredFolder = paste(gseString,  "_untared", sep = "")
celFilesPath = paste("~/a/data", gseUntaredFolder, sep = "/")
setwd(celFilesPath)

library(affy)
library(sva)
pd = read.AnnotatedDataFrame(filename="pheno.txt")
afyData = affy::ReadAffy(phenoData=pd, sampleNames=pd$SampleAccessionNumber)
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
                                         outdir = "QC_affy",
                                         force = TRUE,
                                         do.logtransform = FALSE,
                                         intgroup = "Target")


