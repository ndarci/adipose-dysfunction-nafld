library(ggplot2)
library(ggpubr)

# import data
mt = read.table("../../data/aligned/pass2/mtreadproportions.txt", header = T)

# compute proportions
mt['prop'] = mt$mtreads / mt$totreads
colnames(mt)[1] = "Sample_ID"

# merge with sample/condition information
cond = read.csv("../../data/condition_per_sample_clean.csv", header = T)
names(cond)[1] = 'Sample_ID'
mt = merge(mt, cond, by = "Sample_ID")

# generate bar plot
for (protname in c('CCDC80', 'SOD3')) {
	mtplot = ggbarplot(mt[mt$protein_name == protname,], 
			x = "protein_amt", 
			y = "prop", 
			add = "mean_sd", 
			add_params = list(group = "protein_amt"),
			position = position_dodge(0.8),
			fill = "lightgrey") +
		xlab("Protein amount") + ylab("Proportion MT reads") +
		coord_cartesian(ylim=c(0, 0.25))
	ggsave(paste0("../../figs/mt_read_proportion_bar_", protname, ".png"))	
}

