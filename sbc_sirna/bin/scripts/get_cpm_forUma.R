library(edgeR)

# import cleaned data
cond = read.table("../../data/condition_per_sample_clean.txt", sep = '\t', header = T)
counts = read.table("../../data/knockdown_counts_clean.txt", sep = '\t', header = T)

# read counts and groups into a DGElist
y = DGEList(counts = counts, group = cond$Timepoint_Condition)

# filter for expressed genes
thresh = 3/43 # agreed upon threshold
keepgenes = rowSums(y$counts > 0) >= thresh * length(colnames(y$counts))
table(keepgenes)
y = y[keepgenes, , keep.lib.sizes = F]

# compute CPMs
y = calcNormFactors(y)
write.table(cpm(y), "../../data/cpm_for_uma_allSamples_7pctOver0.txt", quote = F, col.names = T, row.names = T, sep = '\t')





