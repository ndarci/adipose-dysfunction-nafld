library(GSEABase)
library(TissueEnrich)
library(dplyr)
library(ggplot2)
library(SummarizedExperiment)

setwd("~/src/adiDys/bin/scripts")

## for teGeneRetrieval function...
# read in median TPMs from GTEx
expData <- read.table(file = "../../data/deGenesMedianTPM_GTEx.txt")
se <- SummarizedExperiment(assays = SimpleList(as.matrix(expData)), rowData = row.names(expData), colData = colnames(expData))
# run TissueEnrich
output <- teGeneRetrieval(se)
# clean up output
teGenes <- data.frame(assay(output))
teGenes <- teGenes[teGenes$Group %in% c("Tissue-Enhanced", "Group-Enriched"),]
write.table(teGenes, file = "../../data/tissueEnrichedDEGenes.txt", quote = FALSE)
# filter for SBCs
sbc_sym <- read.table("../../data/sbc_IDs.txt")
teGenes_sbc <- teGenes[teGenes$Gene %in% sbc_sym$V1, ]
write.table(teGenes_sbc, file = "../../data/tissueEnrichedDEGenes_sbc.txt", quote = FALSE)

## for teEnrichment function...
# # import background list of all expressed genes
# bg <- scan("../../data/allExprGenes.txt", character())
# gsBG <- GeneSet(geneIds = bg, organism = "Homo Sapiens", geneIdType=ENSEMBLIdentifier())
# import serum biomarker candidates into a GeneSet object
sbc <- scan("../../data/sbc_IDs_ensembl.txt", character())
gsSBC <- GeneSet(geneIds = sbc, organism = "Homo Sapiens", geneIdType=ENSEMBLIdentifier())
# import DE genes into a GeneSet object
de <- scan("../../data/allDEgenes.txt", character())
gsDE <- GeneSet(geneIds = de, organism = "Homo Sapiens", geneIdType=ENSEMBLIdentifier())
# make GeneSet objects for individual genes
gsCD300LG <- GeneSet(geneIds = "ENSG00000161649", organism = "Homo Sapiens", geneIdType=ENSEMBLIdentifier())
gsCPM <- GeneSet(geneIds = "ENSG00000135678", organism = "Homo Sapiens", geneIdType=ENSEMBLIdentifier())
gsCCDC80 <- GeneSet(geneIds = "ENSG00000091986", organism = "Homo Sapiens", geneIdType=ENSEMBLIdentifier())

runTE <- function(gs, title) {
  # run TissueEnrich
  # out <- teEnrichment(inputGenes = gs, backgroundGenes = gsBG)
  out <- teEnrichment(inputGenes = gs)
  # isolate output
  seEnrichmentOutput <- out[[1]]
  eo <- setNames(data.frame(assay(seEnrichmentOutput), row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
  eo$Tissue <- row.names(eo)
  # get table of SE results
  setable <- data.frame(assay(seEnrichmentOutput))
  setable["p-value"] <- 10 ^ -setable["Log10PValue"]
  # plot results
  p <- ggplot(eo,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue, label = Tissue.Specific.Genes,fill = Tissue)) +
    geom_bar(stat = 'identity') +
    labs(x='', y = '-LOG10(P-Value)') +
    theme_bw() +
    theme(legend.position='none') +
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.grid.major= element_blank(),panel.grid.minor = element_blank()) + 
    ggtitle(title)
  retme <- list(tab = setable, plot = p)
  return(retme)
}

teSBC <- runTE(gsSBC, "Tissue enrichment for serum biomarker candidates")
teDE <- runTE(gsDE, "Tissue enrichment for all DE genes")
teCD300LG <- runTE(gsCD300LG, "Tissue enrichment for CD300LG")
teCPM <- runTE(gsCPM, "Tissue Enrichment for CPM")
teCCDC80 <- runTE(gsCCDC80, "Tissue Enrichment for CCDC80")



