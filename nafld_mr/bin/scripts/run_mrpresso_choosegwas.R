library(MRPRESSO)

# import prepped MR data
mr = read.table('../../data/mr_input_table_TwoSampleMR.txt')

# Run MR-PRESSO global method
mrpresso = mr_presso(BetaOutcome = "beta.outcome", 
		BetaExposure = "beta.exposure", 
		SdOutcome = "se.exposure", 
		SdExposure = "se.outcome", 
		OUTLIERtest = TRUE, 
		DISTORTIONtest = TRUE, 
		data = mr, 
		NbDistribution = 1000,  
		SignifThreshold = 0.05)
print(mrpresso)

write.table(mrpresso$`Main MR results`, '../../data/mrpresso_results_choosegwas.txt')
