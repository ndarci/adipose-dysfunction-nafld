tissue = commandArgs(trailingOnly = T)[1]

library(WGCNA)

# set global data.frame options
options(stringsAsFactors = FALSE)

# Read in the expression dataset
lnames = load(file = paste0("../../data/", tissue, "_clean_counts_cov.RData"))
datExpr0 = datExpr
#datExpr0 = read.table(paste0("../../data/", tissue, "_logcpm_covcorrect_INT_topGenes.txt"))

# do some QC checks on the data

# check for missingness
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
# if there was missingness, remove the bad samples
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# hierarchically cluster the samples to find outliers
sampleTree = hclust(dist(datExpr0), method = "average")
# par(cex = 0.6)
# par(mar = c(0,4,2,0))

# plot the clusters
png(paste0("../../fig/qc_clustering_outliers_", tissue, ".png"), height = 7, width = 7, units = 'in', res = 150)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# remove the one outlier we found by cutting the tree
# Plot a line to show the cut
abline(h = 15, col = "red")
dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)

if (tissue == "adipose") {
  # no outliers detected
  keepSamples = (clust==0)
  datExpr = datExpr0[keepSamples, ]
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  } else { # tissue is liver
  keepSamples = (clust==0)
  datExpr = datExpr0[keepSamples, ]
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
}

# make sure expression and phenotype data include the same samples
datExpr = datExpr[rownames(datExpr) %in% rownames(datTraits),]

# make a heatmap relating sample clustering to phenotypes
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
png(paste0("../../fig/qc_clustering_heatmap_", tissue, ".png"), height = 7, width = 7, units = 'in', res = 150)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()




