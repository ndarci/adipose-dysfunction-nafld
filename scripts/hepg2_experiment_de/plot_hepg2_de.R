library(edgeR)
library(tidyverse)
library(gridExtra)
library(gtable)
library(grid)

upcolor = 'dodgerblue2'
downcolor = 'indianred2'
nonsigcolor = 'grey'
controlcolor = 'grey33'

dotplottreatmentcolor = 'deeppink'

theme_set(theme_bw())

amt_order = c('0ng', '20ng')
# prot_order = c('CCDC80')
prot_order = c('CCDC80', 'SOD3')

# read in DE results
read_in_de <- function(nominal = F) {
	de = data.frame()
	for(protname in prot_order) {
		covcorrect = F 
		covariate = ""
		fn = paste0("../../data/topTable_", 
					protname, 
					rep(paste0("_covcorrected_", covariate), covcorrect), 
					rep("_nominal", nominal),
					".txt")
		newdata = read.table(fn, header = T, sep = '\t')
		if (nrow(newdata) > 0) {
			newdata["protein_name"] = protname
			newdata["covcorrect"] = covcorrect
			newdata["covariate"] = covariate
			de = rbind(de, newdata)
		}
	}
	return(de)
}

de = read_in_de(nominal = T)

for (protname in prot_order) {

	# get results from this protein treatment
	de.sub = de %>% filter(protein_name == protname)

	# ##############
	# # generate a volcano plot for each experiment
	# ##############

    # # identify upregulated, downregulated, non-significant
    # de.sub = de.sub %>% 
    # 	mutate(threshold = ifelse(adj.P.Val>0.05, 'ns', 
    # 						ifelse(logFC>0, 'up', 'down')),
    # 			gene_symbol = ifelse(adj.P.Val<0.05, gene_name, ''))

    # de.sig = de.sub %>% filter(adj.P.Val < 0.05) %>%
    # 		mutate(xnudge = 0,
    # 				ynudge = 0)
    # mypos = 'top'

	# assign(paste0('volcano.', protname), ggplot(de.sub, aes(x = logFC,
    #                                     y = -log10(adj.P.Val),
    #                                     label = gene_symbol)) +
    #     geom_point(aes(color = threshold)) +
    #     scale_color_manual(breaks = c('up', 'down', 'ns'),
	# 		            labels = c('Upregulated   ',
    #                             'Downregulated   ',
    #                             'Non-sig. (adjP>0.05)       '),
	# 		            values = c(upcolor, downcolor, nonsigcolor)) +
    #     geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    #     annotate("text", x = 3.6, y = 1.25, label = "adjP=0.05") +
    #     theme(legend.position = mypos,
    #                     legend.title = element_blank(),
    #                     panel.grid.minor = element_blank(),
    #                     text = element_text(size = 14)) +
    #     xlab(paste0("logFC between protein treatment and control")) +
    #     ylab("-log10(Adjusted p-value)") +
    #     ggrepel::geom_text_repel(data = de.sig,
    #                         nudge_x = de.sig$xnudge,
    #                         nudge_y = de.sig$ynudge,
    #                         show.legend = F))
    #     ggsave(paste0("../../figs/volcanoplot_", protname, ".png"), get(paste0('volcano.', protname)), dpi = 300)

	##############
	# generate a logFC bar plot for each experiment
	##############

	# throw away nominal genes
	de.sub = de.sub %>% filter(adj.P.Val < 0.05)

	# sort by logFC
	fc_gene_order = de.sub %>% arrange(desc(logFC)) %>% .$gene_name
	de.sub['gene_name'] = factor(de.sub$gene_name, level = rev(fc_gene_order))

	# define color column based on DE significance and direction
	de.sub['direc_sig'] = as.numeric(de.sub$logFC > 0)
	# special color for non-significant
	de.sub[de.sub$adj.P.Val > 0.05,'direc_sig'] <- 2
	de.sub['direc_sig'] = factor(de.sub$direc_sig, level = c('0', '1'))

	# plot logFC of every DE gene
	assign(paste0('logFCbar.', protname), ggplot(de.sub, aes(x = logFC, 
						y = gene_name,
						fill = direc_sig)) + 
					geom_bar(stat = 'identity') +
					scale_fill_manual(name = '',
						breaks = c('2', '1', '0'),
						labels = c('Non-sig. (adjP>0.05)    ',
									paste0('Up in ', protname, ' treatment  '), 
									paste0('Down in ', protname, ' treatment')),
						values = c(nonsigcolor, upcolor, downcolor),
						drop = F) +
					ylab('') + xlab('logFC between protein treatment and control') +
					xlim(-2, 2.3) + 
					guides(fill = guide_legend(nrow = 2)) +
					theme(panel.grid.minor = element_blank(),
						legend.position = 'top',
						text = element_text(size = 14),
						strip.background = element_rect(fill = 'white')))
	ggsave(paste0('../../figs/logFC_barplot_', protname, '.png'), get(paste0('logFCbar.', protname)), width = 9, height = 8, units = 'in', dpi = 300)

	##############
	# Next part: generate dot and line plots showing the change 
	# in expression with protein treatment
	##############

	# import cleaned expression data
	# crop to only samples used for this protein [CCDC80 or SOD3]
	cond = read.csv("../../data/condition_per_sample_clean.csv") %>%
			filter(protein_name == protname)
	counts = read.table("../../data/hepg2_counts_clean.txt", sep = '\t', header = T) %>%
			select(cond$sample)
	gene_info = read.table("../../data/gene_info.txt", sep = '\t', header = T, )
	gene_info['gene_id'] = rownames(gene_info)

	# import counts and conditions into a DGElist
	y = DGEList(counts = counts, group = cond$protein_amt)

	# remove lowly expressed genes
	MINCOUNT = 10
	keepgenes = filterByExpr(y, min.count = MINCOUNT, group = cond$protein_amt)
	y = y[keepgenes, , keep.lib.sizes = F]

	# compute CPM
	y = calcNormFactors(y)
	cpm = cpm(y) %>% t(.) %>% data.frame(.)

	# restrict to only the genes we care about
	cpm = cpm[,de.sub$gene_ID]
	names(cpm) = de.sub$gene_name

	# merge condition and expression together
	cpm['sample'] = rownames(cpm)
	ov = merge(cond, cpm, by = 'sample')

	# plot mean +-sd expression of control and protein treatment for all DE genes
	ovsub = ov %>% filter(protein_name == protname)

	# standardize the CPM by the mean of each gene to put them on the same axis
	standardize_col <- function(col) {
		newcol = (col - mean(col)) / (sd(col))
		return(newcol)
	}
	ovsub = ovsub %>% mutate_at(names(cpm)[1:length(names(cpm))-1],
				standardize_col)

	# format data for plotting
	ovsub = ovsub %>% pivot_longer(cols = de.sub$gene_name,
									names_to = 'exprgene',
									values_to = 'cpm') %>%
				group_by(protein_amt, exprgene) %>%
				summarise(meanCPM = mean(cpm),
							barlo = mean(cpm)-sd(cpm),
							barhi = mean(cpm)+sd(cpm))

	ovsub = merge(ovsub, de[,c('gene_name', 'logFC', 'P.Value', 'adj.P.Val')], 
		by.x = 'exprgene', 
		by.y = 'gene_name', all.x = T)

	# define color for KD vs scramble
	ovsub['colorcode'] = ovsub$protein_amt == '0ng'
	ovsub['colorcode'] = factor(as.numeric(ovsub$colorcode))

	# sort genes the same way they were in the barplot
	ovsub['exprgene'] = factor(ovsub$exprgene, level = rev(fc_gene_order))

	# put p-values for each DE test in a plottable format
	anno_df = ovsub %>% group_by(exprgene) %>%
					summarise(adj.P.Val = adj.P.Val[!(is.na(adj.P.Val))],
								P.Value = P.Value[!(is.na(P.Value))],
								y.text = mean(c(max(barlo), min(barhi)))) %>%
					unique()

	# set significance annotation text
	anno_df['sig'] = ''
	# anno_df[anno_df$adj.P.Val < 0.1,'sig'] <- '.'
	anno_df[anno_df$adj.P.Val < 0.05,'sig'] <- '*'
	anno_df[anno_df$adj.P.Val < 0.01,'sig'] <- '**'
	anno_df[anno_df$adj.P.Val < 0.001,'sig'] <- '***'

	# generate the plot itself
	assign(paste0('dotplot.', protname), 
					ggplot(ovsub, 
					aes(y = exprgene,
						x = meanCPM,
					color = colorcode)) + 
		scale_color_manual(breaks = c('0', '1'),
						labels = c(paste0(protname, ' treatment   '), 'Control'),
						values = c(dotplottreatmentcolor, controlcolor)) +
		geom_errorbar(aes(xmin = barlo, xmax = barhi), width = 0, alpha = 0.5) +
		geom_point() +
		xlab('Standardized gene expression') +
		xlim(-1.75, 1.75) + 
		ylab('') +
		labs(color = '') +
		guides(color = guide_legend(nrow = 2)) +
		geom_text(data = anno_df, nudge_y = -0.025, 
								aes(y = exprgene,
										x = y.text,
										label = sig,
										group = NULL,
										color = NULL),
									show.legend = F) +
		theme(panel.grid.minor = element_blank(),
			legend.position = 'top',
			text = element_text(size = 14),
			strip.background = element_blank(),
			strip.text = element_blank()))
	ggsave(paste0('../../figs/expression_dotplot_', protname, '.png'), get(paste0('dotplot.', protname)), width = 5.5, height = 12, units = 'in', dpi = 300)
}


# arrange into one huge multipanel with the plot areas lined up

# convert to grobs
c.grob.fc = ggplotGrob(logFCbar.CCDC80)
c.grob.dot = ggplotGrob(dotplot.CCDC80)
s.grob.fc = ggplotGrob(logFCbar.SOD3)
s.grob.dot = ggplotGrob(dotplot.SOD3)

# lay the plots out in a 2x2 grid
layout = rbind(c(1, 2),
				c(3, 4))
both.multi = grid.arrange(c.grob.fc, c.grob.dot, 
							s.grob.fc, s.grob.dot,
				layout_matrix = layout, widths = c(7.5, 6.5),
				heights = c(8.15, 2.85))
ggsave('../../figs/multipanel_both_hepg2.png', both.multi,
	width = 14, height = 12, units = 'in', dpi = 800)










