library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(cowplot)

# import previously generated heatmap plots
lnames = load('../../fig/gg_module_trait_corr_heatmap_adipose.RData')
a.heat = ggcorr
lnames = load('../../fig/gg_module_trait_corr_heatmap_liver.RData')
l.heat = ggcorr
lnames = load('../../fig/gg_crosstissue_ME_correlation_heatmap.RData')
x.heat = ggcorr

# get the legend by itself to put on the side
corr.legend = cowplot::get_legend(x.heat)

# generate a blank panel for the layout
blank = grid.rect(gp = gpar(col = 'white'))

# make the plots into grobs so they can be aligned
a.heat.grob = ggplotGrob(a.heat + theme(axis.text.y = element_blank(),
									legend.position = 'none') +
								scale_y_discrete(limits = a.module_order))
l.heat.grob = ggplotGrob(l.heat + theme(axis.text.x = element_blank(),
									legend.position = 'none') +
								scale_x_discrete(limits = l.module_order))
x.heat.grob = ggplotGrob(x.heat + theme(legend.position = 'none'))

# align left side of liver plot with left side of xtissue plot
l.left = l.heat.grob$widths
x.left = x.heat.grob$widths

max.left = unit.pmax(l.left, x.left)
l.heat.grob$widths = max.left
x.heat.grob$widths = max.left

# align bottom of adipose plot with bottom of xtissue plot
a.bottom = a.heat.grob$heights
x.bottom = x.heat.grob$heights

max.bottom = unit.pmax(a.bottom, x.bottom)
a.heat.grob$heights = max.bottom
x.heat.grob$heights = max.bottom

# define layout matrix of the figure
layout = rbind(c(3, 1, 5),
				c(4, 2, 5))

# arrange all plots into a figure
multi = grid.arrange(blank, 
					a.heat.grob, 
					l.heat.grob, 
					x.heat.grob, 
					corr.legend,
					layout_matrix = layout,
					widths = c(20, 6, 2),
					heights = c(1, 6))
ggsave('../../fig/multipanel_wgnca_heatmap.png', multi, 
	width = 15, height = 18, units = 'in', dpi = 800)











