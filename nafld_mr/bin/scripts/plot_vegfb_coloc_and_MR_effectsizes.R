library(tidyverse)
library(ggrepel)
library(gridExtra)
library(gtable)
library(grid)

##################
# first thing: generate a locuszoom plot for VEGFB colocalized region
##################

# read in SNP overlap data between TG gwas hits and vegfb cis-eQTLs
ov = read.table('../../data/snp_overlap/TG_noNAFLD_ENSG00000173511_SNPoverlap.txt')
ov = ov %>% select(rsID_gwas, chr_eqtl, pos_eqtl, p_eqtl, p_gwas)

# read in local LD information
ld = read.table('../../data/ld_matrix/ld_matrix_ENSG00000173511.ld')
# import SNPs in the LD matrix (overlap of KOBS, UKB, METSIM)
snplist = read.table(paste0('../../data/ld_matrix/ld_matrix_ENSG00000173511.snplist'), header = F)$V1
rssplit = strsplit(snplist, ':', fixed = T)
snplist_rsonly = sapply(rssplit, '[[', 1)
# give the ld matrix its row/col names
rownames(ld) = snplist_rsonly
names(ld) = snplist_rsonly

coolsnp = 'rs2845885'

coolsnp.gwas = 'rs56271783'

# get r2 with the SNP we care about
ld = ld %>% select(coolsnp.gwas) %>% 
		.^2 %>% 
		rename(r2_ivsnp = coolsnp.gwas) %>%
		mutate(rsID_gwas = rownames(ld))

# add to the big table
ov = ov %>% left_join(ld, by = 'rsID_gwas') %>% na.omit(.)

# set up for discrete color scale
colbreaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
colors = c('blue4','skyblue','darkgreen','orange','red')
ov$r2color = as.character(
						cut(ov$r2_ivsnp, 
							breaks = colbreaks, 
							include.lowest = T)
						)
colororder = c('[0,0.2]', '(0.2,0.4]', '(0.4,0.6]', '(0.6,0.8]', '(0.8,1]')

# set up for labeled point
ov = ov %>% mutate(snplabel = ifelse(rsID_gwas == coolsnp, rsID_gwas, ''),
					snplabel.gwas = ifelse(rsID_gwas == coolsnp.gwas, rsID_gwas, ''))

theme_set(theme_bw())

# generate eqtl locuszoom plot
textsize = 14
zoom_eqtl = ggplot(ov, aes(x = pos_eqtl/1e6, 
							y = -log10(p_eqtl), 
							color = r2color)) + 
	geom_point() +
	geom_point(data = ov %>% filter(rsID_gwas == coolsnp.gwas), 
				color = 'purple', shape = 'diamond', size = 3) + 
	geom_text_repel(aes(label = snplabel), color = 'black', 
					nudge_x = -0.35, nudge_y = -0.5,
					segment.linetype = 'dotted') + 
	geom_text_repel(aes(label = snplabel.gwas), color = 'black', 
				nudge_x = 0.5, nudge_y = -0.25,
				segment.linetype = 'dotted') + 
	scale_color_manual(values = colors, breaks = colororder) + 
	labs(x = element_blank(),
		y = expression('VEGFB '*italic('cis')*'-eQTL '-'log'[10]*'(P)')) +
	theme(legend.position = 'none',
		text = element_text(size = textsize))
ggsave('../../figs/custom_locuszoom_eqtl.png', zoom_eqtl)

# generate gwas locuszoom plot
zoom_gwas = ggplot(ov, aes(x = pos_eqtl/1e6, 
							y = -log10(p_gwas), 
							color = r2color)) + 
	geom_point() +
	geom_point(data = ov %>% filter(rsID_gwas == coolsnp.gwas), 
			color = 'purple', shape = 'diamond', size = 3) + 
	geom_text_repel(aes(label = snplabel), color = 'black', 
					nudge_x = -0.4, nudge_y = 1,
					segment.linetype = 'dotted') + 
	geom_text_repel(aes(label = snplabel.gwas), color = 'black', 
			nudge_x = 0.5, nudge_y = -0.5,
			segment.linetype = 'dotted') + 
	scale_color_manual(values = colors, breaks = colororder) + 
	labs(x = 'Chromosome 11 (Mb)', 
		y = expression('Triglycerides GWAS '-'log'[10]*'(P)')) +
	theme(legend.position = 'none',
		text = element_text(size = textsize))
ggsave('../../figs/custom_locuszoom_gwas.png', zoom_gwas)

zoom_both = ggplot(ov, aes(x = -log10(p_gwas),
							y = -log10(p_eqtl), 
							color = r2color)) +
	geom_point() + 
	geom_point(data = ov %>% filter(rsID_gwas == coolsnp.gwas), 
			color = 'purple', shape = 'diamond', size = 3) + 
	geom_text_repel(aes(label = snplabel), color = 'black', 
					nudge_x = -15, nudge_y = -0.25,
					segment.linetype = 'dotted') + 
	geom_text_repel(aes(label = snplabel.gwas), color = 'black', 
			nudge_x = -2, nudge_y = -1,
			segment.linetype = 'dotted') + 
	scale_color_manual(values = colors, breaks = colororder) + 
	guides(color = guide_colorsteps(barwidth = 1,
									show.limits = T,
									title = expression('r'^2),
									title.hjust = 0.25,
									frame.colour = 'black',
									frame.linewidth = 1)) +
	labs(x = expression('Triglycerides GWAS '-'log'[10]*'(P)'),
		y = expression('VEGFB '*italic('cis')*'-eQTL '-'log'[10]*'(P)')) +
	annotate(geom = 'text', x = 38, y = 11.75, label = paste0('LD=', round(ld[ld$rsID_gwas == coolsnp,'r2_ivsnp'], 2))) +
	theme(legend.position = c(0.85, 0.175),
			legend.box.background = element_rect(color = 'black'),
			text = element_text(size = textsize))
ggsave('../../figs/custom_locuszoom_both.png', zoom_both)


# plot multipanel locuszoom plot
zoom.multi = grid.arrange(zoom_both, zoom_eqtl, zoom_gwas, 
				layout_matrix = rbind(c(1, 2),
									c(1, 3)))
ggsave('../../figs/multipanel_locuszoom_vegfb.png', zoom.multi, width = 7.7, height = 7, units = 'in', dpi = 800)


##################
# second thing: generate a scatterplot showing MR effect sizes for TG and NAFLD
##################


# import prepped MR data
mr = read.table('../../data/mr_input_table_TwoSampleMR.txt')

# import MR-PRESSO results
mrpresso = read.table('../../data/mrpresso_results.txt')

# import other 3 method results
othermr = read.table('../../data/other_mr_methods_results.txt')

# plot the MR results
# trim to what we need, and add columns for se
coolsnp = 'rs2845885'
mr = mr %>% 
		select(SNP, beta.exposure, beta.outcome, se.exposure, se.outcome) %>%
		mutate(
			se.lo.exposure = beta.exposure - se.exposure,
			se.hi.exposure = beta.exposure + se.exposure,
			se.lo.outcome = beta.outcome - se.outcome,
			se.hi.outcome = beta.outcome + se.outcome,
			snplabel = ifelse(SNP == coolsnp, SNP, ''))

ndig = 2
mrlabel = paste0('MR-PRESSO: beta = ', round(mrpresso[1,'Causal.Estimate'], ndig), ';',
				' P = ', formatC(mrpresso[1,'P.value'], format = 'e', digits = ndig), '\n',
				'MR-PRESSO horizontal pleiotropy test: P = NS')

mrplot = ggplot(mr, aes(x = beta.exposure, y = beta.outcome)) +
		geom_abline(intercept = 0, 
					slope = mrpresso[1,'Causal.Estimate'],
					color = 'red', alpha = 1) +
		geom_errorbarh(aes(xmin = se.lo.exposure, xmax = se.hi.exposure), height = 0, color = 'gray22') +
		geom_errorbar(aes(ymin = se.lo.outcome, ymax = se.hi.outcome), width = 0, color = 'gray22') +
		geom_text_repel(aes(label = snplabel), color = 'black',
						nudge_x = -0.02, nudge_y = 0.002,
						segment.linetype = 'dotted') +
		geom_point() + 
		xlab('Beta for triglycerides') + ylab('Beta for imputed NAFLD') +
		annotate(geom = 'label', x = 0.013, y = -0.0035, label = mrlabel, hjust = 0) +
		theme(text = element_text(size = 14))
ggsave('../../figs/mrpresso_scatterplot.png', mrplot)


# put both the plots together into a full size figure
zoom.both.grob = ggplotGrob(zoom_both + 
				theme(axis.title.y = element_text(vjust = -9, size = 14)))
zoom.gwas.grob = ggplotGrob(zoom_gwas)
zoom.eqtl.grob = ggplotGrob(zoom_eqtl)
mr.grob = ggplotGrob(mrplot)

left.zoom = zoom.both.grob$widths[1:4]
left.mr = mr.grob$widths[1:4]

max.left = unit.pmax(left.zoom, left.mr)
zoom.both.grob$widths[1:4] = max.left
mr.grob$widths[1:4] = max.left


layout.full = rbind(c(1, 2),
				c(1, 3),
				c(4, 4),
				c(4, 4))

coloc.mr.plot = grid.arrange(zoom.both.grob, zoom.eqtl.grob, 
							zoom.gwas.grob, mr.grob,
							layout_matrix = layout.full)
ggsave('../../figs/coloc_mrpresso_multi.png', coloc.mr.plot, 
		width = 8, height = 14, units = 'in', dpi = 800)











