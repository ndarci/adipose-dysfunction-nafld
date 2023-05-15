library(tidyverse)
library(gridExtra)
library(cowplot)
library(gtable)
library(grid)
library(stringr)
theme_set(theme_bw())

# iterate over all the wgcna modules we care about
for (thing in list(c('liver', 'cyan'),
					c('liver', 'darkmagenta'),
					c('liver', 'royalblue'),
					c('liver', 'saddlebrown'),
					c('liver', 'tan'),
					c('liver', 'violet'),
					c('adipose', 'lightyellow'),
					c('adipose', 'cyan'))) {
	tissue = thing[[1]]
	module = thing[[2]]

	# import pathway enrichment results from webgestalt
	filepath = paste0('../../data/enrichment_results_wg_result_', tissue, '_', module, '.txt')
	enrich = read.table(filepath, sep = '\t', header = T)

	# fix weird double autophagy
	if (tissue == 'liver' & module == 'royalblue') {
		enrich[10, 'description'] = 'Autophagy_1'
	}
	if (tissue == 'liver' & module == 'cyan') {
		enrich[10, 'description'] = 'Glycosphingolipid biosynthesis_1'
	}

	if (length(unique(enrich$description)) != length(enrich$description)) {
		print(paste('fuck', tissue, module))
	}

	# choose columns we want
	enrich = enrich %>% 
			select(description, size, overlap, expect, 
				enrichmentRatio, pValue, FDR) %>%
			mutate(signif = FDR < 0.05) %>%
			arrange(enrichmentRatio) %>%
			mutate(description = factor(description, level = unique(description)))

	# make a bar plot with this information
	sigcolor = 'black'
	nonsigcolor = 'grey'
	assign(paste0('plot.', tissue, '.', module), 
					ggplot(enrich, aes(x = enrichmentRatio, 
										y = description,
										fill = signif)) +
					geom_bar(stat = 'identity') +
					scale_fill_manual(name = '',
						breaks = c('TRUE', 'FALSE'),
						labels = c('Significant (FDR<0.05)     ',
									'Non-significant (FDR>0.05)'),
						values = c(sigcolor, nonsigcolor)) +
					xlab('Enrichment Ratio') +
					ggtitle(paste0(str_to_title(tissue), ' ', module)) +
					theme(axis.title.y = element_blank(),
						legend.position = 'none',
						plot.title = element_text(hjust = 0.5),
						text = element_text(size = 14))
		)
	assign(paste0('grob.', tissue, '.', module),
				ggplotGrob(get(paste0('plot.', tissue, '.', module))))
	# ggsave(paste0('../../fig/enrichment_barplot_', tissue, '_', module, '.png'), 
	# 		enrichbarplot, width = 10, height = 4, units = 'in', dpi = 800)	
}

# grab the legend from a plot with some sig and some not sig
plot.liver.violet = plot.liver.violet + theme(legend.position = 'top')
masterlegend = cowplot::get_legend(plot.liver.violet)
plot.liver.violet = plot.liver.violet + theme(legend.position = 'none')

# ggsave('../../fig/test.png', plot.adipose.lightyellow)

# arrange all the bar plots into one big multipanel figure
layout = rbind(c(10, 9, 9, 11),
				c(10, 1, 2, 11),
				c(10, 3, 4, 11),
				c(10, 5, 6, 11),
				c(10, 7, 8, 11))

# set right col width of plot areas
right.maxwidth = unit.pmax(grob.adipose.cyan$widths,
							grob.liver.tan$widths,
							grob.liver.darkmagenta$widths,
							grob.liver.cyan$widths)
grob.adipose.cyan$widths = right.maxwidth
grob.liver.tan$widths = right.maxwidth
grob.liver.darkmagenta$widths = right.maxwidth
grob.liver.cyan$widths = right.maxwidth
# set left col width of plot areas
# do same as right ones to keep sizes consistent

# left.maxwidth = unit.pmax(grob.adipose.lightyellow$widths,
# 							grob.liver.saddlebrown$widths,
# 							grob.liver.violet$widths,
# 							grob.liver.royalblue$widths)

grob.adipose.lightyellow$widths = right.maxwidth
grob.liver.saddlebrown$widths = right.maxwidth
grob.liver.violet$widths = right.maxwidth
grob.liver.royalblue$widths = right.maxwidth


multi = grid.arrange(grob.adipose.lightyellow, grob.adipose.cyan,
					grob.liver.saddlebrown, grob.liver.tan,
					grob.liver.violet, grob.liver.darkmagenta,
					grob.liver.royalblue, grob.liver.cyan,
					masterlegend,
					layout_matrix = layout,
					heights = c(2, 10, 10, 10, 10),
					widths = c(1, 10, 10, 1))

ggsave('../../fig/multipanel_enrichment_barplot.png', multi,
		height = 18, width = 16, units = 'in', dpi = 800)




