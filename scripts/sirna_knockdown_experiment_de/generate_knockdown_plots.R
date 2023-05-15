library(edgeR)
library(tidyverse)
library(gridExtra)
library(gtable)
library(grid)

time_order = c("Baseline", "24h", "4D", "7D")
sbc_order = c("CCDC80", "SOD3")

kdcolor = 'darkgoldenrod1'
upcolor = 'dodgerblue2'
downcolor = 'indianred2'
nonsigcolor = 'grey'
controlcolor = 'grey33'

lineplotkdcolor = 'deeppink'

theme_set(theme_bw())

# plot logFC of DE genes at each timepoint
read_in_de <- function(nominal = F) {
	# import DE results
	de = data.frame()
	for (tp in time_order) {
		for (sbc in sbc_order) {
			# for (covariate in c(paste0("SV1_", sbc), "unique_map_pct", "")) {
			# 	if (covariate == "") { 
					covcorrect = F 
					covariate = ""
				# } else { covcorrect = T }
				fn = paste0("../../data/topTable_", 
							tp, "_", 
							sbc, 
							rep(paste0("_covcorrected_", covariate), covcorrect), 
							rep("_nominal", nominal),
							".txt")
				newdata = read.table(fn, header = T, sep = '\t')
				if (nrow(newdata) > 0) {
					newdata["timepoint"] = tp
					newdata["knockdown"] = sbc
					newdata["covcorrect"] = covcorrect
					newdata["covariate"] = covariate
					de = rbind(de, newdata)
				# }
			}
		}
	}
	return(de)
}

de = read_in_de(nominal = T)
de['timepoint'] = factor(de$timepoint, level = time_order)

for (kd in sbc_order) {
	# get genes from this knockdown
	de.sub = de[de$knockdown == kd,]

	# define set of genes for the plot: those that are significant at any timepoint
	siganytime = de.sub[de.sub$adj.P.Val < 0.05,]$gene_name %>% unique(.)
	de.sub = de.sub[de.sub$gene_name %in% siganytime,] 

	# sort genes by logFC on important day
	if (kd == 'CCDC80') {sortday = '7D'} else {sortday = '4D'}
	fc_gene_order = de.sub %>% filter(timepoint == sortday) %>% arrange(desc(logFC)) %>% .$gene_name
	de.sub['gene_name'] = factor(de.sub$gene_name, level = rev(fc_gene_order))

	# define color column based on DE significance and direction
	de.sub['direc_sig'] = as.numeric(de.sub$logFC > 0)
	# special color for non-significant
	de.sub[de.sub$adj.P.Val > 0.05,'direc_sig'] <- 2
	# special color for the knockdown itself
	de.sub[de.sub$gene_name == kd & de.sub$P.Value < 0.05, 'direc_sig'] <- 3
	de.sub['direc_sig'] = factor(de.sub$direc_sig)

	# plot logFC per timepoint
	assign(paste0('logFCbar.', kd), ggplot(de.sub, aes(x = logFC, 
						y = gene_name,
						fill = direc_sig)) + 
					geom_bar(stat = 'identity') +
					facet_grid(.~timepoint) +
					scale_fill_manual(name = '',
						breaks = c('3', '2', '1', '0'),
						labels = c(paste0(kd, ' KD  '),
									'Non-sig. (adjP>0.05)    ',
									paste0('Up in ', kd, ' KD  '), 
									paste0('Down in ', kd, ' KD')),
						values = c(kdcolor, nonsigcolor, upcolor, downcolor)) +
					ylab('') + xlab('logFC between KD and scramble') +
					guides(fill = guide_legend(nrow = 2)) +
					theme(panel.grid.minor = element_blank(),
						legend.position = 'top',
						text = element_text(size = 14),
						strip.background = element_rect(fill = 'white')))
	ggsave(paste0('../../figs/logFC_barplot_', kd, '.png'), get(paste0('logFCbar.', kd)), width = 9, height = 8, units = 'in', dpi = 300)
}

# import cleaned expression data
cond = read.table("../../data/condition_per_sample_clean.txt", sep = '\t', header = T)
descols = c('Sample_ID', 'Sample_Name', 'Timepoint', 'Condition', 'Timepoint_Condition', 'Replicate')
cond = cond[,descols]
counts = read.table("../../data/knockdown_counts_clean.txt", sep = '\t', header = T)
gene_info = read.table("../../data/gene_info.txt", sep = '\t', header = T, )
gene_info['gene_id'] = rownames(gene_info)

# define genes we want to look at in the line plots
c.targets = c('SREBF1', 'SCD', 'LPL', 'AHR')
s.targets = c('NHLRC3', 'LEP', 'G0S2', 'RBL1')
coolgenes = data.frame(knockdown = c(rep('CCDC80', length(c.targets)), rep('SOD3', length(s.targets)), rep(NA, length(sbc_order))),
						gene_name = c(c.targets, s.targets, sbc_order))
coolgenes = merge(coolgenes, gene_info[,c('gene_name', 'gene_id')], by = 'gene_name', all.x = T)

# remove weird RBL1 double ID
coolgenes = coolgenes %>% filter(gene_id != 'ENSG00000269846.1')

# import counts and conditions into a DGElist
y = DGEList(counts = counts, group = cond$Condition)

# remove lowly expressed genes
MINCOUNT = 10
keepgenes = filterByExpr(y, min.count = MINCOUNT, group = cond$Timepoint_Condition)
y = y[keepgenes, , keep.lib.sizes = F]

# compute CPM
y = calcNormFactors(y)
cpm = cpm(y) %>% t(.) %>% data.frame(.)

# restrict to only the genes we care about
cpm = cpm[,coolgenes$gene_id]
names(cpm) = coolgenes$gene_name
coolgenes = na.omit(coolgenes)

# merge condition and expression together
cpm['Sample_ID'] = rownames(cpm)
ov = merge(cond, cpm, by = 'Sample_ID')

# plot mean +-sd expression of control and knockdown for key target genes over time
for (kd in c('CCDC80', 'SOD3')) {
	# select desired target genes
	target = coolgenes[coolgenes$knockdown == kd,]$gene_name
	if (kd == 'CCDC80') {
		target_order = c.targets 
		coolday = '7D'
	} else {
		target_order = s.targets
		coolday = '4D'
	}
	# format data for plotting
	ovsub = ov[ov$Condition == 'Controli' | ov$Condition == kd, c(descols, kd, target)]
	ovsub = ovsub %>% pivot_longer(cols = c(kd, target),
									names_to = 'exprgene',
									values_to = 'cpm') %>%
				group_by(Timepoint, Condition, exprgene) %>%
				summarise(meanCPM = mean(cpm),
							barlo = mean(cpm)-sd(cpm),
							barhi = mean(cpm)+sd(cpm))
	ovsub['Timepoint'] = factor(ovsub$Timepoint, level = time_order)
	ovsub['exprgene'] = factor(ovsub$exprgene, level = c(kd, target))
	ovsub['Condition'] = factor(ovsub$Condition, level = c('Controli', kd))
	ovsub['exprgene'] = factor(ovsub$exprgene, level = c(kd, target_order))

	# add DE stats for color/annotation
	ovsub = merge(ovsub, de[,c('timepoint', 'knockdown', 'gene_name', 
								'logFC', 'P.Value', 'adj.P.Val')], 
			by.x = c('Timepoint', 'Condition', 'exprgene'), 
			by.y = c('timepoint', 'knockdown', 'gene_name'), all.x = T)

	# # define color for each line corresponding to direction of that gene on most important day
	# gene_color = ovsub[ovsub$Timepoint == coolday,]
	# gene_color['colorcode'] = as.numeric(gene_color$logFC > 0)
	# # set knocked down gene to a different color
	# gene_color[gene_color$exprgene == kd & gene_color$Condition == kd, 'colorcode'] <- 3
	# # set non-significant to a different color
	# gene_color[is.na(gene_color$colorcode), 'colorcode'] <- 2
	# # merge color info back onto data table
	# ovsub = merge(ovsub, gene_color[,c('Condition', 'exprgene', 'colorcode')], 
	# 				by = c('Condition', 'exprgene'))

	# define color for KD vs scramble
	ovsub['colorcode'] = ovsub$Condition == 'Controli'
	ovsub['colorcode'] = factor(as.numeric(ovsub$colorcode))

	# put p-values for each DE test in a plottable format
	anno_df = ovsub %>% group_by(exprgene, Timepoint) %>%
					summarise(adj.P.Val = adj.P.Val[!(is.na(adj.P.Val))],
								P.Value = P.Value[!(is.na(P.Value))],
								y.text = mean(meanCPM))

	# set significance annotation text for non-KD genes
	anno_df['sig'] = ''
	# anno_df[anno_df$adj.P.Val < 0.1,'sig'] <- '.'
	anno_df[anno_df$adj.P.Val < 0.05,'sig'] <- '*'
	anno_df[anno_df$adj.P.Val < 0.01,'sig'] <- '**'
	anno_df[anno_df$adj.P.Val < 0.001,'sig'] <- '***'

	# set special annotations for KD genes
	anno_df[anno_df$exprgene == kd & anno_df$P.Value < 0.05, 'sig'] <- '+'
	anno_df[anno_df$exprgene == kd & anno_df$P.Value < 0.01, 'sig'] <- '++'
	anno_df[anno_df$exprgene == kd & anno_df$P.Value < 0.001, 'sig'] <- '+++'

	# nudge pesky stars
	if(kd == 'CCDC80') {
		anno_df[anno_df$exprgene == 'LPL' &
				anno_df$Timepoint == '4D','y.text'] <- 800
	} else if (kd == 'SOD3') {
		anno_df[anno_df$exprgene == 'SOD3' &
				anno_df$Timepoint == 'Baseline', 'y.text'] <- 35
		anno_df[anno_df$exprgene == 'LEP' &
				anno_df$Timepoint == '4D', 'y.text'] <- 100
		anno_df[anno_df$exprgene == 'G0S2' &
				anno_df$Timepoint == '24h', 'y.text'] <- 1100
		anno_df[anno_df$exprgene == 'RBL1' &
				anno_df$Timepoint == '4D', 'y.text'] <- 3		
	}

	assign(paste0('kdplot.', kd), ggplot(ovsub %>% arrange(Condition), aes(x = Timepoint, 
						y = meanCPM,
						color = colorcode,
						group = Condition)) + 
			# scale_color_manual(breaks = c('3', '2', '1', '0'),
			# 					labels = c(paste0(kd, ' KD  '),
			# 								'Scramble',
			# 								paste0('Up in ', kd, ' KD  '), 
			# 								paste0('Down in ', kd, ' KD  ')),
			# 					values = c(kdcolor, controlcolor, upcolor, downcolor)) +
			scale_color_manual(breaks = c('0', '1'),
							labels = c(paste0(kd, ' KD   '), 'Scramble'),
							values = c(lineplotkdcolor, controlcolor)) +
			facet_grid(exprgene~., scales = 'free_y') +
			geom_errorbar(aes(ymin = barlo, ymax = barhi), width = 0, alpha = 0.5) +
			geom_line() +
			geom_point() +
			ylab('CPM') +
			labs(color = '') +
			geom_text(data = anno_df, aes(x = Timepoint,
											y = y.text,
											label = sig,
											group = NULL,
											color = NULL),
										show.legend = F) +
			guides(color = guide_legend(nrow = 2)) +
			theme(panel.grid.minor = element_blank(),
				legend.position = 'top',
				text = element_text(size = 14),
				strip.background = element_rect(fill = 'white')))
	ggsave(paste0('../../figs/knockdown_lines_', kd, '.png'), get(paste0('kdplot.', kd)), width = 5.5, height = 8, units = 'in', dpi = 300)


}

# # arrange these plots in a multipanel figure
# c.multi = grid.arrange(logFCbar.CCDC80, kdplot.CCDC80, 
# 					nrow = 1, widths = c(3, 2))
# ggsave('../../figs/multipanel_CCDC80_knockdown.png', c.multi, 
# 	width = 14, height = 8, units = 'in', dpi = 800)

# s.multi = grid.arrange(logFCbar.SOD3, kdplot.SOD3,
# 					nrow = 1, widths = c(3, 2))
# ggsave('../../figs/multipanel_SOD3_knockdown.png', s.multi, 
# 	width = 14, height = 8.5, units = 'in', dpi = 800)

# both.multi = grid.arrange(logFCbar.CCDC80, kdplot.CCDC80,
# 							logFCbar.SOD3, kdplot.SOD3,
# 							nrow = 2, widths = c(3, 2), heights = c(8, 8.5))
# ggsave('../../figs/multipanel_both_knockdown.png', both.multi, 
# 	width = 14, height = 16.5, units = 'in', dpi = 800)


# arrange into one huge multipanel with the plot areas lined up

# convert to grobs
c.grob.fc = ggplotGrob(logFCbar.CCDC80)
s.grob.fc = ggplotGrob(logFCbar.SOD3)
c.grob.kd = ggplotGrob(kdplot.CCDC80)
s.grob.kd = ggplotGrob(kdplot.SOD3)

# set left col width of plot areas
left.maxwidth = unit.pmax(c.grob.fc$widths, s.grob.fc$widths)
c.grob.fc$widths = left.maxwidth
s.grob.fc$widths = left.maxwidth
# same for right col
right.maxwidth = unit.pmax(c.grob.kd$widths, s.grob.kd$widths)
c.grob.kd$widths = right.maxwidth
s.grob.kd$widths = right.maxwidth

# # same for top and bottom rows
# top.maxheight = unit.pmax(c.grob.fc$heights, c.grob.kd$heights)
# c.grob.fc$heights = top.maxheight
# c.grob.kd$heights = top.maxheight

# lay the plots out in a 2x2 grid
layout = rbind(c(1, 2),
				c(3, 4))
both.multi = grid.arrange(c.grob.fc, c.grob.kd, s.grob.fc, s.grob.kd, 
						widths = c(3, 2), heights = c(8, 9),
						layout_matrix = layout)
ggsave('../../figs/multipanel_both_knockdown.png', both.multi, 
	width = 14, height = 17, units = 'in', dpi = 800)



