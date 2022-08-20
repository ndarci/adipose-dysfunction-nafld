tissue = commandArgs(trailingOnly = T)[1]

library(WGCNA)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Load the normalized expression data
lnames = load(file = paste0("../../data/", tissue, "_clean_counts_cov.RData"))

# pick the soft threshold value used to define network connectivity
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results
cex1 = 0.9
png(paste0("../../fig/soft_thresh_by_r2_", tissue, ".png"), height = 7, width = 7, units = 'in', res = 150)
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
dev.off()

# Mean connectivity as a function of the soft-thresholding power
png(paste0("../../fig/soft_thresh_by_mean_connectivity_", tissue, ".png"), height = 7, width = 7, units = 'in', res = 150)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


if (tissue == "adipose") {
     # we choose this threshold b/c it's the lowest threshold where the  
     # topology model fit is high and mean connectivity is leveling off
     softPower = 7
     # decide on a ME cut height threshold of 0.25 == 0.75 correlation
     MEDissThres = 0.10
} else {
     # tissue is liver
     softPower = 10
     MEDissThres = 0.25
}

# calculate the adjacency matrix
adjacency = adjacency(datExpr, power = softPower)

# convert adj matrix into topological overlap matrix
# this is a better measure of gene relatedness than pure adjacency
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# compute hierarch. clustering on the dissimilarity TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
png(paste0("../../fig/initial_clustering_", tissue, ".png"), height = 7, width = 7, units = 'in', res = 150)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# identify modules by applying a dynamic tree cut to this dendrogram
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
png(paste0("../../fig/initial_clustering_withcolors_", tissue, ".png"), height = 7, width = 7, units = 'in', res = 150)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# merge modules with highly similar expression profiles, using module eigengenes
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result


png(paste0("../../fig/ME_clustering_", tissue, ".png"), height = 7, width = 7, units = 'in', res = 150)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()


# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# plot the new merged modules
png(paste0("../../fig/merged_clustering_withcolors_", tissue, ".png"), height = 7, width = 7, units = 'in', res = 150)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# save the final merged data for downstream
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = paste0("../../data/", tissue, "_constructed_network.RData"))
