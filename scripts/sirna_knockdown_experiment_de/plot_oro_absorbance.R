library(tidyverse)
library(rstatix)
theme_set(theme_bw())

controlcolor = 'grey33'
lineplotkdcolor = 'deeppink'

# read in the ORO absorbance data used to measure adipocyte lipidation
oro = read.csv('../../data/oro_absorbance_table.csv', row.names = 1)

# clean it all up
time_order = c('Baseline', '24h', '4D', '7D')
oro = oro %>% 
		# swap rows for columns
		t() %>%
		data.frame() %>%
		# add clean timepoint column
		mutate(timepoint = c(rep('Baseline', 4),
							rep('24h', 4),
							rep('4D', 4),
							rep('7D', 4))) %>%
		# # add control column, all 1s bc relative to this reading
		# mutate(Control = 1) %>%
		# pivot longer for plot formatting
		pivot_longer(cols = c('Ci', 'CCDC80', 'SOD3'),
					names_to = 'condition',
					values_to = 'absorbance') %>%
		# make sure timepoints go in the right order
		mutate(timepoint = factor(timepoint, levels = time_order))

# compute t-tests comparing scramble to knockdown at each timepoint
df.annot = data.frame()
for (t in time_order) {
	addme = oro %>% 
		filter(timepoint == t) %>% 
		pairwise_t_test(absorbance ~ condition) %>%
		filter(group1 == 'Ci' | group2 == 'Ci') %>%
		mutate(timepoint = t)
	df.annot = rbind(df.annot, addme)
}
df.annot = df.annot %>% 
			select(timepoint, group1, group2, p) %>%
			mutate(kd = ifelse(group1 != 'Ci', group1, group2),
				sig = ifelse(p<0.001, '***', 
						ifelse(p<0.01, '**',
							ifelse(p<0.05, '*', ''))),
				y.text = ifelse(kd == 'CCDC80', 0.75, 0.825))

# copy the control data so we can have it for CCDC80 and SOD3 plots
oro.noci = oro %>% filter(condition != 'Ci')
oro.ci = oro %>% filter(condition == 'Ci')
oro.full = rbind(oro.noci %>% mutate(kd = condition),
				oro.ci %>% mutate(kd = 'CCDC80'),
				oro.ci %>% mutate(kd = 'SOD3'))

# summarise into mean and sd per timepoint/condition
oro.sum = oro.full %>%
		group_by(timepoint, condition, kd) %>%
		summarise(meanABS = mean(absorbance),
				barlo = mean(absorbance) - sd(absorbance),
				barhi = mean(absorbance) + sd(absorbance)) %>%
		mutate(colorcode = condition != 'Ci',
			condition = recode(condition, Ci = 'Scramble'))

# plot the absorbance over time for scramble, ccdc80, and sod3 against control
oroplot = ggplot(oro.sum, aes(x = timepoint, 
								y = meanABS,
								group = colorcode,
								color = colorcode)) +
			scale_color_manual(breaks = c(TRUE, FALSE),
							labels = c('Knockdown', 'Scramble'),
							values = c(lineplotkdcolor, controlcolor)) +
			facet_grid(kd~., scales = 'free_y') +
			geom_errorbar(aes(ymin = barlo, ymax = barhi), width = 0, alpha = 0.5) +
			geom_line() + 
			geom_point() +
			xlab('Timepoint') + 
			ylab('Absorbance relative to non-scramble control') +
			labs(color = '') +
			guides(color = guide_legend(nrow = 2)) +
			theme(panel.grid.minor = element_blank(),
				legend.position = 'top',
				text = element_text(size = 14),
				strip.background = element_rect(fill = 'white')) +
			geom_text(data = df.annot, aes(x = timepoint,
											y = y.text,
											label = sig,
											group = NULL,
											color = NULL),
										show.legend = F)

ggsave('../../figs/oro_absorbance_plot.png', oroplot)









