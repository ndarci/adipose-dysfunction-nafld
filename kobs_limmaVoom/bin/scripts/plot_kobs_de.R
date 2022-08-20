library(tidyverse)
library(gridExtra)
library(cowplot)
library(gtable)
library(grid)

theme_set(theme_bw())

upcolor = 'dodgerblue2'
downcolor = 'indianred2'
nonsigcolor = 'grey'

hist_order = c('steatosis', 'fibrosis', 'diagnosis')
hist_map = list(steatosis = 'steatosis',
				fibrosis = 'fibrosis',
				diagnosis = 'NASH')

# read in kobs DE results
read_in_de <- function(tissue = 'adipose', nominal = F) {
	de = data.frame()
	for (pheno in hist_order) {
		fn = paste0("../../data/topTable_",
					pheno, "_",
					tissue,
					rep("_nominal", nominal),
					".txt")
		newdata = read.table(fn, header = T)
		if (nrow(newdata) > 0) {
			newdata["phenotype"] = pheno
			newdata["gene_ID"] = rownames(newdata)
			rownames(newdata) <- NULL
			de = rbind(de, newdata)
		}
	}
	return(de)
}
de = read_in_de(nominal = T)

# import SBC names
sbctable <- read.table("../../data/metTable_subset_serumBiomarkers.txt", sep = '\t', header = TRUE)

# create volcano plots for each phenotype
for (p in hist_order) {
	# select results for this phenotype
	desub = de[de$phenotype == p,]

	# merge SBC names onto table
	desub = merge(desub, unique(sbctable[sbctable$phenotype == p, c('gene_ID', 'gene_symbol')]), by = 'gene_ID', all.x = T)

	# identify upregulated, downregulated, non-significant
	desub = desub %>% mutate(threshold = ifelse(adj.P.Val>0.05, 'ns', ifelse(logFC>0, 'up', 'down')))
	desbc = desub[!is.na(desub$gene_symbol),]
	print(desbc)

	# define label nudges
	if(p == 'steatosis') {
		#				CCDC80, TIMP3, SOD3, COL6A1, COL6A2, SRFP2, GPX3
		desbc['xnudge'] = c(0.3, -0.6, -0.9, 1, 1, 0.8, -0.6)
		desbc['ynudge'] = c(2, 2, 0.1, -1, 0.1, 0.3, -1)
		mypos = 'none'
	} else if (p == 'fibrosis') {
		#				CCDC80, SOD3, SFRP2
		desbc['xnudge'] = c(0.1, -0.3, 0.2)
		desbc['ynudge'] = c(1.3, 1.1, 1.2)
		mypos = 'none'
	} else {
		# p is nash
		# 				CCDC80, TIMP3, SOD3, MGP, SRFP2, CD300LG, VEGFB
		desbc['xnudge'] = c(0.1, -1, -0.4, 0.6, 0.8, -0.8, -0.8)
		desbc['ynudge'] = c(1.5, 0.8, 2, -0.9, 0.8, -0.3, 0.3)
		mypos = 'top'
	}

	assign(paste0('volcano.', p), ggplot(desub, aes(x = logFC, 
								y = -log10(adj.P.Val),
								label = gene_symbol)) +
				geom_point(aes(color = threshold)) +
				scale_color_manual(breaks = c('up', 'down', 'ns'),
									labels = c('Upregulated   ',
												'Downregulated   ',
												'Non-sig. (adjP>0.05)       '),
									values = c(upcolor, downcolor, nonsigcolor)) +
				geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
				annotate("text", x = 1.3, y = 1.2, label = "adjP=0.05") +
				theme(legend.position = mypos,
						legend.title = element_blank(),
						panel.grid.minor = element_blank(),
						text = element_text(size = 14)) +
				xlab(paste0("Adipose logFC between ", hist_map[p], " and healthy liver")) + 
				ylab("-log10(Adjusted p-value)") +
				ggrepel::geom_text_repel(data = desbc,
										nudge_x = desbc$xnudge,
										nudge_y = desbc$ynudge,
										show.legend = F))
	ggsave(paste0("../../figs/volcanoplot_", p, ".png"), get(paste0('volcano.', p)), dpi = 300)
}
volc.legend = cowplot::get_legend(volcano.diagnosis)

# plot inverted bars of logFC per phenotype for all SBCs
# get subset of SBC results for all DE tests
defc = de[de$gene_ID %in% sbctable$gene_ID,]
defc = merge(defc, sbctable[,c('gene_ID', 'gene_symbol')], by = 'gene_ID')

# define color column based on DE significance and direction
defc['direc_sig'] = as.numeric(defc$logFC > 0)
defc[defc$adj.P.Val > 0.05,'direc_sig'] <- 2
defc['direc_sig'] = factor(defc$direc_sig)

# recode to make the plot labels prettier
defc['phenotype'] = factor(defc$phenotype, level = hist_order) %>%
							recode(.,
							steatosis = 'Steatosis',
							fibrosis = 'Fibrosis',
							diagnosis = 'NASH')

# sort according to logFC
defc = defc %>% arrange(desc(logFC))
defc['gene_symbol'] = factor(defc$gene_symbol, level = rev(unique(defc$gene_symbol)))


# plot logFC per timepoint
logFCbar = ggplot(defc, 
				aes(x = logFC, 
					y = gene_symbol,
					fill = direc_sig)) + 
				geom_bar(stat = 'identity') +
				facet_grid(.~phenotype) +
				scale_fill_manual(name = '',
					breaks = c('1', '0', '2'),
					labels = c(paste0('Upregulated'), paste0('Downregulated'), 'Non-significant (adjP>0.05)'),
					values = c(upcolor, downcolor, nonsigcolor)) +
				ylab('') + xlab('Adipose logFC between sick and healthy liver') +
				theme(panel.grid.minor = element_blank(),
					legend.position = 'top',
					text = element_text(size = 14),
					strip.background = element_rect(fill = 'white'))
logFC.legend = cowplot::get_legend(logFCbar)
ggsave('../../figs/logFC_barplot_allphenotypes.png', logFCbar, dpi = 300)

# multi = grid.arrange(volcano.steatosis, volcano.fibrosis, 
# 					volcano.diagnosis + theme(legend.position = 'none'), 
# 					volc.legend,
# 					logFCbar + theme(legend.position = 'none'), logFC.legend, 
# 					layout_matrix = rbind(c(1, 2, 3, 4),
# 										c(5, 5, 5, 6)),
# 					widths = c(5, 5, 5, 1),
# 					heights = c(2, 1))

# convert to grobs so plot areas match up
nash.grob = ggplotGrob(volcano.diagnosis + theme(axis.title.y = element_text(vjust = -12, size = 14)))
logFC.grob = ggplotGrob(logFCbar + theme(legend.position = 'none'))
# align left facet with left side of big plot, and same with right side
left.nash = nash.grob$widths[1:4]
left.logFC = logFC.grob$widths[1:4]

max.left = unit.pmax(left.nash, left.logFC)
nash.grob$widths[1:4] = max.left
logFC.grob$widths[1:4] = max.left

multi = grid.arrange(nash.grob,
					logFC.grob,
					nrow = 2)
ggsave('../../figs/multipanel_KOBS_DE_nashvolcano.png', multi,
	width = 7.5, height = 13, units = 'in', dpi = 800)


supp.multi = grid.arrange(
							volcano.steatosis + theme(legend.position = 'top'),
							volcano.fibrosis + theme(legend.position = 'none'),
							nrow = 2)
ggsave('../../figs/multipanel_volcano_steatosis_fibrosis.png', supp.multi,
	width = 7, height = 13, units = 'in', dpi = 800)




