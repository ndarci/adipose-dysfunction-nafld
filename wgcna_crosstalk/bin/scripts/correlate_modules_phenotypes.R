tissue = commandArgs(trailingOnly = T)[1]

library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# Load the expression and trait data saved in the first part
lnames = load(file = paste0("../../data/", tissue, "_clean_counts_cov.RData"))
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = paste0("../../data/", tissue, "_constructed_network.RData"))
lnames

# quantify assoc of module eigengenes with phenotypes
# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# remember number of genes per module
# get names (colors) of the modules
modNames = substring(names(MEs), 3)
me_count = data.frame(table(moduleColors))
colnames(me_count) = c("module", "count")
rownames(me_count) = me_count$module
me_count = me_count[modNames,]
modNames_count = paste0(modNames, " (", me_count$count, ")")

# compute correlations of modules and traits
# first, remove covariates (we know they have no correlation)
covlist = c("Age", "isMale", "RIN", "UNIQ_MAP_PERCENT", 
          "PCT_INTRONIC_BASES", "MEDIAN_3PRIME_BIAS")
datTraits = datTraits[,!(names(datTraits) %in% c(covlist, "Steatosisgradebaseline", "Fibrosisstagebaseline", "Diagnosisbaseline"))]
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# prep for bonferroni correction on the p-values
nModules = dim(moduleTraitPvalue)[1]
nPheno = dim(moduleTraitPvalue)[2]
adj_p_thresh = 0.05 / (nModules * nPheno)

# # display correlations and their p-values
# textMatrix = moduleTraitPvalue < adj_p_thresh
# textMatrix[textMatrix == T] <- ""
# textMatrix[textMatrix == F] <- "X"
# # Display the correlation values within a heatmap plot
# png(paste0("../../fig/module_trait_corr_heatmap_", tissue, ".png"), height = 9, width = 5, units = 'in', res = 150)
# par(mar = c(7, 8, 1, 1))
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels = names(datTraits),
#                yLabels = paste0("ME", modNames_count),
#                ySymbols = names(MEs),
#                colorLabels = FALSE,
#                colors = blueWhiteRed(50),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.4,
#                cex.lab = 0.75,
#                zlim = c(-1,1))
# dev.off()


# make a custom heatmap with ggplot
library(reshape2)
library(tidyverse)

# calculate first pc for ordering
pc1 = prcomp(moduleTraitCor)$x[,1] %>% 
        data.frame(.) %>%
        rename('pc1' = 1) %>%
        mutate('module' = rownames(.))

# set up to have gene counts with each module
df_mod_count = data.frame(module = rownames(moduleTraitCor),
                            module_count = modNames_count) %>%
                mutate(module_count = paste0(toupper(substr(tissue, 0, 1)), 
                                            substr(tissue, 2, nchar(tissue)),
                                            ' ', 
                                            module_count)) %>%
                inner_join(pc1, by = 'module') %>%
                arrange(pc1)

# get data into correct format
            # pivot r values into tidy format
cor.melt = melt(moduleTraitCor) %>%
            # merge with p values
            inner_join(melt(moduleTraitPvalue), by = c('Var1', 'Var2')) %>%
            # clean up column names
            rename(module = 'Var1',
                    trait = 'Var2',
                    r = 'value.x',
                    p = 'value.y') %>%
            # add name plus count column
            inner_join(df_mod_count, by = 'module') %>%
            # add column for significance annotation
            mutate(nonsig = as.character(p > adj_p_thresh) %>%
                    recode('TRUE' = 'X', 'FALSE' = '')) %>%
            # make phenotype names pretty
            mutate(trait = factor(trait, 
                    levels = c('dm_1.yes', 'BMI', 'cholesterolmed_preoper', 
                                'Fasting_glucose', 'Triglycerides',
                                'steatosis_grp', 'fibrosis_grp', 'diagnosis_grp'))) %>%
            mutate(trait = recode(trait, 
                    'dm_1.yes' = 'T2D',
                    'steatosis_grp' = 'Steatosis',
                    'fibrosis_grp' = 'Fibrosis',
                    'diagnosis_grp' = 'NASH',
                    'cholesterolmed_preoper' = 'Cholesterol',
                    'Fasting_glucose' = 'Fasting glucose'))

# generate the plot
if (tissue == 'adipose') {
    ggcorr = ggplot(cor.melt, aes(x = trait, y = module_count, fill = r)) +
        geom_tile() + 
        geom_text(aes(label = nonsig)) +
        # coord_fixed() +
        scale_fill_gradient2(low = 'darkorchid2', 
                            mid = 'white', 
                            high = 'forestgreen',
                            limits = c(-1, 1)) +
        guides(fill = guide_colorbar(title = expression(italic('R')),
                                    barwidth = 0.5, 
                                    barheight = 15,
                                    label.position = 'right',
                                    title.hjust = 0,
                                    ticks = F)) + 
        scale_y_discrete(limits = df_mod_count$module_count, position = 'left') +
        theme(panel.background = element_blank(),
                legend.position = 'right',
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                text = element_text(size = 14))
} else {
    ggcorr = ggplot(cor.melt, aes(x = module_count, y = trait, fill = r)) +
        geom_tile() + 
        geom_text(aes(label = nonsig)) +
        # coord_fixed() +
        scale_fill_gradient2(low = 'darkorchid2', 
                            mid = 'white', 
                            high = 'forestgreen',
                            limits = c(-1, 1)) +
        guides(fill = guide_colorbar(title = expression(italic('R')),
                                    barwidth = 0.5, 
                                    barheight = 15,
                                    label.position = 'right',
                                    title.hjust = 0,
                                    ticks = F)) + 
        scale_x_discrete(limits = df_mod_count$module_count, position = 'top') +
        theme(panel.background = element_blank(),
                legend.position = 'right',
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.x = element_text(angle = 45, vjust = -1, hjust = -1),
                text = element_text(size = 14))
}

save(ggcorr, file = paste0('../../fig/gg_module_trait_corr_heatmap_', tissue, '.RData'))
ggsave(paste0('../../fig/gg_module_trait_corr_heatmap_', tissue, '.png'), 
    ggcorr, width = 8, height = 12, units = 'in', dpi = 800)






# # compute gene significance and module membership
# # Define variable coolPheno containing some column of datTrait
# coolPheno = as.data.frame(datTraits$steatosis_grp)
# names(coolPheno) = "steatosis_grp"

# # compute membership
# geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
# names(geneModuleMembership) = paste("MM", modNames, sep="")
# # compute p values of membership
# MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
# names(MMPvalue) = paste("p.MM", modNames, sep="")

# # compute significance
# geneTraitSignificance = as.data.frame(cor(datExpr, coolPheno, use = "p"))
# names(geneTraitSignificance) = paste("GS.", names(coolPheno), sep="")
# GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
# names(GSPvalue) = paste("p.GS.", names(coolPheno), sep="")

# # visualize GS x MM in the module most assoc with coolPheno
# module = substring(names(which.max(moduleTraitCor[,"steatosis_grp"])), 3)

# column = match(module, modNames)
# moduleGenes = moduleColors==module
# # sizeGrWindow(7, 7)
# png(paste0("../../fig/membership_vs_signif_", tissue, "_", names(coolPheno), ".png"), height = 7, width = 7, units = 'in', res = 150)
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = paste0("Gene significance for ", names(coolPheno)),
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# dev.off()

# # match gene IDs to gene names
# annot = read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/gencodeV26_annotations_formatted_genomewide.txt", header = TRUE, sep = '\t')
# mygenes = data.frame(names(datExpr))
# names(mygenes) = "gene_id"
# annot = merge(mygenes, unique(annot[,c("gene_id", "gene_name")]), all.x = T)
# annot = annot[match(names(datExpr), annot$gene_id),]

# # package all these results nicely into a dataframe
# # Create the starting data frame
# geneInfo0 = data.frame(gene_ID = names(datExpr),
#                        gene_symbol = annot$gene_name,
#                        moduleColor = moduleColors,
#                        geneTraitSignificance,
#                        GSPvalue)
# # Order modules by their significance for coolPheno
# modOrder = order(-abs(cor(MEs, coolPheno, use = "p")))
# # Add module membership information in the chosen order
# for (mod in 1:ncol(geneModuleMembership))
# {
#   oldNames = names(geneInfo0)
#   geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
#                          MMPvalue[, modOrder[mod]]);
#   names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
#                        paste("p.MM.", modNames[modOrder[mod]], sep=""))
# }
# # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.steatosis_grp))
# geneInfo = geneInfo0[geneOrder, ]

# write.csv(geneInfo0, file = paste0("../../data/", tissue, "_geneInfo.csv"))
