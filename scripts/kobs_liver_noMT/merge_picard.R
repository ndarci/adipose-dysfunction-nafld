# extract tech factors from picard output files for all liver samples

# get sample IDs
samples = read.table("sampleIDs.txt", header = F)$V1

# set up an empty dataframe to paste each row onto
picout = data.frame()

# iterate through all samples
for (s in samples) {
	# point to picard output file and read its relevant line of data
	picfile = paste("kobs_liver_bam_noMT/", s, "/picard_RNAmetrics_noMT.txt", sep = "")
	pic = read.table(picfile, skip = 6, nrows = 1, header = T, sep = '\t')
	# keep tract of what sample this is
	pic["SAMPLE"] = s
	# merge onto the table with every other sample
	picout = rbind(picout, pic)
}

# write the final table out
write.table(picout, "picardRNAmetrics_merged_kobs_liver_noMT.txt", quote = F, row.names = F, col.names = T, sep = '\t')
