library(ggplot2)
library(ggpubr)

# import data
mt = read.table("../../data/aligned/pass2/mtreadproportions.txt", header = T)

# compute proportions
mt['prop'] = mt$mtreads / mt$totreads
colnames(mt)[1] = "Sample_ID"

# merge with time information
cond = read.table("../../data/condition_per_sample_clean.txt", sep = '\t', header = T)
mt = merge(mt, cond, by = "Sample_ID")

# generate bar plot
mtplot = ggbarplot(mt, 
		x = "Timepoint", 
		y = "prop", 
		add = "mean_sd", 
		add_params = list(group = "Timepoint"),
		position = position_dodge(0.8),
		fill = "lightgrey") +
	xlab("Timepoint") + ylab("Proportion MT reads") +
	coord_cartesian(ylim=c(0, 0.25))
ggsave("../../figs/mt_read_proportion_bar.png")
