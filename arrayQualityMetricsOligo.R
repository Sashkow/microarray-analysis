source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("arrayQualityMetrics")
biocLite("oligo")
biocLite("Biobase")
library(oligo)


gseString = "GSE53291"

gseUntaredFolder = paste(gseString,  "_untared", sep = "")
celFilesPath = paste("~/a/data", gseUntaredFolder, sep = "/")
celFilesPath
setwd(celFilesPath)


pd = Biobase::read.AnnotatedDataFrame(filename = "covdesc")

filePaths = paste(celFilesPath, rownames(pd), sep = "/")
oligoData = oligo::read.celfiles(filenames = filePaths, phenoData = pd)

normalizedData = oligo::rma(oligoData)

arrayQualityMetrics::arrayQualityMetrics(expressionset = normalizedData,
                                         outdir = "QC_oligo",
                                         force = TRUE,
                                         intgroup = "Target")

