# Annotate all CpGs located in the data set according to the parameters below and save as pdf in out output folder.
# ChIPseeker is an R package for annotating ChIP-seq data analysis. 
# It supports annotating ChIP peaks and provides functions to visualize ChIP peaks coverage over chromosomes and profiles of peaks binding to TSS regions. 
# Comparison of ChIP peak profiles and annotation are also supported. Moreover, it supports evaluating significant overlap among ChIP-seq datasets. 
# Currently, ChIPseeker contains 17,000 bed file information from GEO database. 
# These datasets can be downloaded and compare with user’s own data to explore significant overlap datasets for inferring co-regulation or transcription factor complex for further investigation.

# source > https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

# **Note:** Remember to create a "data/" and a "output/" folder in the current directory.

################################################################################

BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("ChIPseeker")

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPseeker)
################################################################################

# Annotate all CpGs located in the data set according to the parameters below and save as pdf in out output folder.

wd <- getwd()
path_to_input= file.path(wd, "data")
path_to_output= file.path(wd, "output")
pattern=".txt" #.narrowPeak
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
samplefiles <- list.files(path = path_to_input, pattern = pattern, full.names = T)
samplefiles <- as.list(samplefiles)
filenames <- gsub(file.path(wd, "data", "*"), "", samplefiles) 
filenames <- gsub("*.narrowPeak", "", filenames)
names(samplefiles) <- filenames
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, verbose = F, tssRegion=c(-1000, 100))
pdf(file = paste0(path_to_output, "/Cluster4_Cdca7_CpGs_with_values_1000-100.pdf"), height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
################################################################################

# Save a table for annotated CpGs as xls or csv

SaveTable <- function(n){
  anno <- data.frame(peakAnnoList[[n]]@anno)
  write.table(anno, file = paste0(wd, "/output/", filenames[n], "_anno.xls"), sep = "\t")}
  #write.table(anno, file = paste0(wd, "/output/", filenames[n], "_anno.csv"), sep = "\t")}

lapply(1:length(samplefiles), SaveTable)

# Calculate the overall frequency in data set
peakAnnoList
################################################################################


