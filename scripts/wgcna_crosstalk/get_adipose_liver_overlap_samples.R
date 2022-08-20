library(magrittr)

# read in count and covariate data
adicountfile = "/u/project/pajukant/nikodm/kobs_limmaVoom/data/rawcounts/BaselineAll_KOBS_geneCount_noMT.txt"
adicounts <- read.table(adicountfile, header = TRUE, row.names = 1)
livcountfile = "/u/project/pajukant/nikodm/kobs_limmaVoom/data_liver/gene_counts.match.txt"
livcounts <- read.table(livcountfile, header = TRUE, row.names = 1)

adicov <- read.csv("/u/project/pajukant/nikodm/kobs_limmaVoom/data/cov.csv", header = TRUE, row.names = 1)
livcov <- read.csv("/u/project/pajukant/nikodm/kobs_limmaVoom/data_liver/cov_liver.csv", header = TRUE, row.names = 2)

# define covariates we care about
covlist = c("Age", "isMale", "RIN", "UNIQ_MAP_PERCENT", "PCT_INTRONIC_BASES", "MEDIAN_3PRIME_BIAS")

# select samples with full data
adicov = na.omit(adicov[,covlist])
livcov = na.omit(livcov[,covlist])

# output list of samples with count and covariate data in both tissues
validsamples = intersect(rownames(adicov), rownames(livcov)) %>%
				intersect(., colnames(adicounts)) %>%
				intersect(., colnames(livcounts))

write.table(validsamples, 
	"../../data/KOBS_adipose_liver_overlap_samples_RNAseq.txt", 
	quote = F, row.names = F, col.names = F)
