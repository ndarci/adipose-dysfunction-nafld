library(data.table)

samples = fread("../../data/sample_fastq.txt")

# define what we mean by polyA
polya = strrep('A', 10)

# keep track of polya counts	
polyacount = c()
polyapct = c()

# iterate over all samples
for (row in seq(1, nrow(samples))) {
	fastqfile = samples[row,]$path
	# read fastq file
	fq = fread(fastqfile, header = F)
	# extract just the reads
	reads = fq[seq(2, nrow(fq), 4),]$V1
	# find reads that start or end with polyA
	startreads = grepl(paste0("^", polya), reads)
	endreads = grepl(paste0(polya, "$"), reads)
	polyareads = !!(startreads + endreads)
	# count polyA reads
	n_polya = sum(polyareads)
	# add to the overall list
	polyacount = c(polyacount, n_polya)
	polyapct = c(polyapct, (n_polya / length(reads)))
}

samples['polyacount'] = polyacount
samples['polyapct'] = polyapct

fwrite(samples[,c("sample", "polyacount", "polyapct")], "../../data/polyAcount_persample.txt")
