setwd("/home/sashko/a/data/")
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("arrayQualityMetrics")
biocLite("simpleaffy")
biocLite("affy")
biocLite("oligo")
gseString = "GSE24129"
gseUntaredFolder = paste(gseString,  "_untared", sep = "")
celFilesPath = paste(getwd(), gseUntaredFolder, sep = "/")
celFilesPath
# --- simpleaffy way ---
# cel files to afffyBatсh
celfiles = simpleaffy::read.affy(path = "GSE24129_untared")

# --- oligo way ---
library(oligo)

# create list of cel files to get raw data from
celFiles <- list.celfiles(celFilesPath, full.names=TRUE)
celFiles

# --- illumnia way ---
library(annotate)
library(beadarray)
library(limma)
library(illuminaMousev2.db)

## The data file from GenomeStudio
dataFile = paste(celFilesPath, "GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt")

## The control file from GenomeStudio
qcFile="/path/to//ControlProbeProfile.txt"

## Reading in the file and creating ExpressionSetIllumina object
eset=readBeadSummaryData(dataFile=dataFile, qcFile=qcFile, 
                         ProbeID="PROBE_ID", controlID="ProbeID", 
                         skip=0, qc.skip=0,
                         annoCols=c("TargetID", "PROBE_ID","SYMBOL"))


# get pheno data 
# from file made by hand according to ncbi site
phenoData <- read.table(paste(celFilesPath, "covdesc", sep = "/"),
                                    header=TRUE,
                                    row.name="Name",
                                    sep="\t")

rawData <- read.celfiles(filenames = celfiles, sampleNames = row.names(phenoData))
pData(rawData) <- phenoData

# before normalization
arrayQualityMetrics::arrayQualityMetrics(expressionset = rawData,
                                         force = TRUE,
                                         do.logtransform = TRUE,
                                         intgroup = "Target")

# --- affy way ---
# celfiles = read.celfiles(p)
# affyBatch to expressionSet
expressionSet = celfiles.rma <- affy::rma(celfiles)

# --- oligo way ---
normalizedRawData2 = backgroundCorrect(rawData)

# after normalization
arrayQualityMetrics::arrayQualityMetrics(expressionset = expressionSet,
                                         force = TRUE,
                                         intgroup = "Target")


design <- model.matrix(~0 + rawData@phenoData@data$Target)

library(limma)
# fit the linear model to the filtered expression set
fit <- lmFit(normalizedRawData2@assayData$exprs, design)

huvec_ebFit <- eBayes(fit)

# return the top 10 results for any given contrast
# coef=1 is huvec_choroid, coef=2 is huvec_retina
topTable(huvec_ebFit, number=10, coef=1)
