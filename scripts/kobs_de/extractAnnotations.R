## extractAnnotations.R
## read the GENCODE v26 primary assembly annotations into a dataframe,
## and extract the information we want to add to the metadata table

#setwd("~/src/adiDys/bin/scripts/")

# read in the raw gtf file
gtf <- read.table("/u/project/pajukant/nikodm/sbc_sirna/data/gencode.v19.annotation.gtf", header = FALSE, sep = '\t')
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# filter for only genes, no transcripts or exons
gtf <- gtf[gtf$feature == "gene", ]

# separate "attribute" column into desired variables
s <- strsplit(gtf$attribute, ";", fixed = TRUE)
m <- Map(function(i, j, k, l) s[[i]][c(j, k, l)],
    i = which(lapply(s, function(x) grep("gene_id", x)) == 1),
    j = unlist(lapply(s, function(x) grep("gene_id", x))),
    k = unlist(lapply(s, function(x) grep("gene_type", x))),
    l = unlist(lapply(s, function(x) grep("gene_name", x))))
library(data.table)
mdf <- transpose(data.frame(m))
colnames(mdf) <- c("gene_id_ver", "gene_type", "gene_name")

mdf$gene_id_ver <- substr(mdf$gene_id_ver, 9, nchar(mdf$gene_id_ver))
mdf$gene_type <- substr(mdf$gene_type, 12, nchar(mdf$gene_type))
mdf$gene_name <- substr(mdf$gene_name, 12, nchar(mdf$gene_name))
sgid <- strsplit(mdf$gene_id_ver, ".", fixed = TRUE)
mdf["gene_id"] <- sapply(sgid, "[[", 1)

# bind extracted variables back onto gtf table
gtf_filt <- cbind(gtf, mdf)
keepcols <- c("gene_id", "gene_id_ver", "seqname", "feature", "start", "end", "strand", "gene_type", "gene_name")
gtf_filt <- gtf_filt[ , keepcols]
gtf_filt <- unique(gtf_filt)

write.table(gtf_filt, file = "../../data/gencodeV19_annotations_formatted_genomewide.txt", quote = FALSE, row.names = FALSE, sep = '\t')

# filter for the genes we care about
myGenes <- read.table("../../data/allDEgenes.txt", header = FALSE)
colnames(myGenes) <- c("gene_id")
gtf_filt <- gtf_filt[gtf_filt$gene_id %in% myGenes$gene_id, ]

# which genes are present in my DE gene table but not in the gtf?
setdiff(myGenes$gene_id, gtf_filt$gene_id)
# what kinds of genes do we have?
table(gtf_filt$gene_type)

# write filtered gtf to tsv
write.table(gtf_filt, file = "../../data/gencodeV19_annotations_filtered.txt", quote = FALSE, row.names = FALSE, sep = '\t')



