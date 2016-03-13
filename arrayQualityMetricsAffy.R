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
affyData = affy::ReadAffy(phenoData=pd,
                          sampleNames=pd$SampleAccessionNumber,
                          cdfname = "nugohs1a520180_hs_entrezg")

length(affyData@assayData$exprs) # probes amnt

exprs <- affy::rma(affyData)
length(exprs@assayData$exprs) # probesets amnt

exprs_frame <- data.frame(exprs)

library(AnnotationDbi)
probesetsID<-rownames(exprs)
length(probesetsID) # amount of annotated probesets in library 17451
probesetsID_EntrezID<-select(nugohs1a520180hsentrezg.db, probesetsID, "ENTREZID")
length(probesetsID_EntrezID$PROBEID) # amount of annotated probesets in library 17451

all <- merge(probesetsID_EntrezID, exprs_frame, by.x=0, by.y=0, all=T)

View(annot)
View(all)


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

arrayQualityMetrics::arrayQualityMetrics(expressionset = exprs,
                                         outdir = "QC_affy_custom_cdf",
                                         force = TRUE,
                                         do.logtransform = FALSE,
                                         intgroup = "Target")

strEndsWith <- function(haystack, needle)
{
  hl <- nchar(haystack)
  nl <- nchar(needle)
  if(nl>hl)
  {
    return(F)
  } else
  {
    return(substr(haystack, hl-nl+1, hl) == needle)
  }
}



