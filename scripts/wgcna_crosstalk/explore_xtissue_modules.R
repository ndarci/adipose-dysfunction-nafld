library(tidyverse)

stripVersion <- function(idlist) {
        strid = strsplit(idlist, ".", fixed = T)
        stripped = sapply(strid, "[[", 1)
        return(stripped)
}

lnames = load('../../data/xtissue_correlation_results.RData')

a.coolmodules = c('lightyellow', 'cyan')
l.coolmodules = c('saddlebrown')

# look further into cross-tissue associated modules
annot = read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/gencodeV26_annotations_formatted_genomewide.txt", sep = '\t', header = T)

a.moduleassign = merge(a.moduleassign, annot[,c('gene_id', 'gene_name')], by = 'gene_id', all.x = T)
l.moduleassign = merge(l.moduleassign, annot[,c('gene_id', 'gene_name')], by = 'gene_id', all.x = T)

# save gene lists for xtissue modules
for (row in seq(1, nrow(sig.pairs))) {
  a.xtissuemodule = substring(sig.pairs$adipose[row], 5)
  l.xtissuemodule = substring(sig.pairs$liver[row], 5)
  a.xtissuegenes = a.moduleassign[a.moduleassign$adipose_module == a.xtissuemodule,]$gene_id
  l.xtissuegenes = l.moduleassign[l.moduleassign$liver_module == l.xtissuemodule,]$gene_id
  write.table(a.xtissuegenes, paste0("../../data/adipose_module_genes_A_", a.xtissuemodule, ".txt"), quote = F, row.names = F, col.names = F)
  write.table(l.xtissuegenes, paste0("../../data/liver_module_genes_L_", l.xtissuemodule, ".txt"), quote = F, row.names = F, col.names = F)
}
# save gene list for statin correlated module
write.table(l.moduleassign[l.moduleassign$liver_module == 'steelblue',]$gene_id, '../../data/liver_module_genes_L_steelblue.txt', quote = F, row.names = F, col.names = F)

# check module membership of cool pathway genes
get_module_mem <- function(genelist, tissuecode, module) {
  if (tissuecode == 'a') tissue = 'adipose' else tissue = 'liver'
  # get ensembl ID
  df = data.frame(gene_name = genelist, 
                  row.names = annot[annot$gene_name %in% genelist,'gene_id'])

  # get module membership for this module
  memberdf = paste0(tissuecode, '.geneModuleMembership')
  df = merge(df, get(memberdf)[paste0('MM_', toupper(tissuecode), '_', module)], by = 0)
  rownames(df) = df$Row.names
  df['Row.names'] = NULL

  # compute membership rank of these genes for all genes and genes in this module 

  # allgenememberships = get(memberdf)[paste0('MM_', toupper(tissuecode), '_', module)]
  # allgenememberships['rank_allgenes'] = rank(allgenememberships[1])
  # allgenememberships['percentile_allgenes'] = allgenememberships$rank_allgenes / nrow(allgenememberships)
  modulegenes = get(paste0(tissuecode, '.moduleassign'))[get(paste0(tissuecode, '.moduleassign'))[paste0(tissue, '_module')] == module,]$gene_id
  thismodulememberships = get(memberdf)[rownames(get(memberdf)) %in% modulegenes,][paste0('MM_', toupper(tissuecode), '_', module)]
  thismodulememberships['rank_thismodule'] = rank(thismodulememberships[1])
  thismodulememberships['percentile_thismodule'] = thismodulememberships$rank_thismodule / nrow(thismodulememberships)
  # df = merge(df, allgenememberships)
  # rownames(df) = df$Row.names
  # df['Row.names'] = NULL
  df = merge(df, thismodulememberships)
  rownames(df) = df$Row.names
  df['Row.names'] = NULL

  return(df)
}
lipolysis = c('LIPE', 'AQP7', 'PLIN1', 'INSR', 'PNPLA2')
lipo_df = get_module_mem(lipolysis, 'a', 'lightyellow')
write.table(lipo_df, '../../data/lipolysis_A_lightyellow_membership.txt', quote = F, row.names = F)
aasynth = c('CPS1', 'CTH', 'ARG1', 'GOT1', 'ACO1', 'ASL', 'ASS1',
            'SDS', 'SDSL', 'MAT1A')
aa_df = get_module_mem(aasynth, 'l', 'saddlebrown')
write.table(aa_df, '../../data/aasynthesis_L_saddlebrown_membership.txt', quote = F, row.names = F)

aly.tf = c('ELF4', 'GTF2E2', 'TCF23', 'IRX6', 'ZBTB7B', 'REPIN1')
write.table(get_module_mem(aly.tf, 'a', 'lightyellow'), '../../data/panther_tf_A_lightyellow_membership.txt', quote = F, row.names = F)
lsb.tf = c('ID2', 'ETV3', 'ZNF281', 'TEAD1', 'ELOA')
write.table(get_module_mem(lsb.tf, 'l', 'saddlebrown'), '../../data/panther_tf_L_saddlebrown_membership.txt', quote = F, row.names = F)

# look at what the shared genes are across lightyellow/saddlebrown
lygenes = a.moduleassign[a.moduleassign$adipose_module == 'lightyellow',]$gene_id
sbgenes = l.moduleassign[l.moduleassign$liver_module == 'saddlebrown',]$gene_id
ly_sb_overlap = data.frame(gene_id = intersect(lygenes, sbgenes))
ly_sb_overlap = merge(ly_sb_overlap, annot[c('gene_name', 'gene_id')], by = 'gene_id')
ly_sb_overlap

# see where my DE genes/SBCs are in important modules
# read in DE gene lists
# these are genes DE in kobs adipose for liver histology, 
# NOT DE in kobs liver, and >10x more expressed in adipose than liver
hist_order = c("steatosis_grp", "fibrosis_grp", "diagnosis_grp")
de = data.frame()
for (hist in hist_order) {
  myGenes <- read.table(paste0("/u/project/pajukant/nikodm/kobs_limmaVoom/data/DE_genes_filter_except_secreted_", hist, ".txt"), header = FALSE)
  myGenes['phenotype'] = hist
  de = rbind(de, myGenes)
}
colnames(de)[1] = 'gene_id'

# test enrichment/overlap of cool modules with DE genes
de = merge(de, a.moduleassign, by = 'gene_id', all.x = T) %>%
      merge(., l.moduleassign, by = 'gene_id', all.x = T) %>%
      merge(., annot[,c("gene_name", "gene_id")], by = "gene_id", all.x = T)

# plot number of DE genes in each adipose/liver module
for (hist in hist_order) {
  a.modcount = data.frame(table(de[de$phenotype == hist,]$adipose_module)) %>% arrange(desc(Freq))
  a.modcount_bar = ggplot(a.modcount, aes(x = factor(Var1, levels = Var1), y = Freq)) +
                    geom_bar(stat = "identity") +
                    xlab("Adipose module") +
                    ylab("Number of adipose-specific DE genes") +
                    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(paste0("../../fig/num_de_genes_permodule_adipose_", hist, ".png"), a.modcount_bar)

  l.modcount = data.frame(table(de[de$phenotype == hist,]$liver_module)) %>% arrange(desc(Freq))
  l.modcount_bar = ggplot(l.modcount, aes(x = factor(Var1, levels = Var1), y = Freq)) +
                    geom_bar(stat = "identity") +
                    xlab("Liver module") +
                    ylab("Number of adipose-specific DE genes") +
                    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(paste0("../../fig/num_de_genes_permodule_liver_", hist, ".png"), l.modcount_bar)
}

# check module membership of SBC genes
sbc = read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/sbc_IDs_ensembl.txt")$V1
de_sub = de[de$gene_id %in% sbc,]
de_sub

# test for de gene enrichment in adipose modules
test_de_enrichment <- function(tissuecode, module) {
  if (tissuecode == 'a') tissue = 'adipose' else tissue = 'liver'
  # count number of genes that are de and in the module
  modulede = de[de[paste0(tissue, '_module')] == module,]$gene_name %>% unique(.)
  nmodulede = length(modulede)
  # count number of genes that are de in the first place
  nde = de$gene_name %>% unique(.) %>% length(.)
  # count number of genes in the whole module
  nmodule = sum(get(paste0(tissuecode, '.moduleassign'))[paste0(tissue, '_module')] == module)
  # get enrichment ratio
  enrich_ratio = (nmodulede / nmodule) / (nde / get(paste0(tissuecode, '.nGenes')))
  # test if there is an enrichment of de genes in this module
  pval_de_enrich = phyper(nmodulede, nde, get(paste0(tissuecode, '.nGenes'))-nde, nmodule, lower.tail = F)

  return(c(nmodulede, nde, nmodule, enrich_ratio, pval_de_enrich, modulede))
}

a.enrich_result = data.frame()
for (m in a.modNames) {
  res = test_de_enrichment('a', m)
  de_overlap_genes = res[6:length(res)]
  if (m %in% a.coolmodules) {
    write.table(de_overlap_genes, paste0('../../data/de_overlap_adipose_', m, '.txt'), quote = F, row.names = F, col.names = F)
  }
  a.enrich_result = rbind(a.enrich_result, data.frame(
                    module = m,
                    n_de_in_module = res[1],
                    n_de = res[2],
                    n_in_module = res[3],
                    enrich_ratio = res[4],
                    pval_de_enrich = res[5]))
}
a.enrich_result['pval_de_enrich'] = as.numeric(a.enrich_result$pval_de_enrich)
write.table(a.enrich_result[a.enrich_result$pval_de_enrich < 0.05,] %>% arrange(pval_de_enrich), '../../data/adipose_de_enrichment.txt', row.names = F)

# import cell type marker data
# adipose
a.marker = read.table('/u/project/pajukant/nikodm/sbc_sirna/data/Cryo.snRNA2019.SAT.n15_markers.celltype_padj0.05_uniq.txt', sep = '\t', header = T)
# add ensembl ID
a.marker['gene_name'] = a.marker$gene
a.marker = merge(a.marker, annot[,c('gene_id', 'gene_name')], by = 'gene_name')

# # liver
# l.marker0 = read.csv('/u/project/pajukant/nikodm/wgcna_crosstalk/data/TableS1.cell_type_markers.csv')
# # select unique cell type marker genes
# l.marker = l.marker0 %>% group_by(Symbol) %>%
#     summarise(Subcell.type = Subcell.type[length(Subcell.type) == 1],
#               Ensembl.ID = Ensembl.ID[length(Subcell.type) == 1],
#               Adjusted.p.value = Adjusted.p.value[length(Subcell.type) == 1]) %>%
#       rename(cluster = Subcell.type,
#           gene_name = Symbol,
#           gene_id = Ensembl.ID,
#           p_val_adj = Adjusted.p.value)
# l.marker['gene_id'] = stripVersion(l.marker$gene_id)

l.marker = read.table('../../data/seur.CellType.markers.unique.txt', header = T)
l.marker['gene_id'] = stripVersion(l.marker$gene)

# run hypergeometric tests to test for enrichment of marker genes
test_ctm_enrichment <- function(tissuecode, module) {
  if (tissuecode == 'a') {
      tissue = 'adipose' 
      celltypes = c('Adipocytes', 'Preadipocytes')
      coolmod = module %in% a.coolmodules
    } else {
      tissue = 'liver'
      celltypes = paste0('Hep-', 1:14)
      coolmod = module %in% l.coolmodules
    }
  # retrieve cell type marker data
  ctm = get(paste0(tissuecode, '.marker'))
  ma = get(paste0(tissuecode, '.moduleassign'))

  # count number of genes in the whole module
  modulegenes = ma[ma[paste0(tissue, '_module')] == module,]$gene_id %>% unique(.)
  nmodule = length(modulegenes)

  modulegenes_sym = ma[ma[paste0(tissue, '_module')] == module,]$gene_name

  thismodule_enrich = data.frame()
  for (c in celltypes) {
    # count number of marker genes for this cell type
    nctm = nrow(na.omit(ctm[ctm$cluster == c,]))
    # count overlap of this cell type marker genes and module genes
    modulectm = intersect(ctm[ctm$cluster == c,]$gene_name, modulegenes_sym) %>% unique(.)
    nmodulectm = length(modulectm)

    # calculate enrichment ratio
    ratio_ctm = (nmodulectm / nmodule) / (nctm / get(paste0(tissuecode, '.nGenes')))
  
    # calculate pvalue of enrichment
    pval_ctm_enrich = phyper(nmodulectm, nctm, get(paste0(tissuecode, '.nGenes'))-nctm, nmodule, lower.tail = F)

    # test enrichment of CTM x DE genes
    ctm_de = intersect(de$gene_name, ctm[ctm$cluster == c,]$gene_name) %>% unique(.)
    # count genes that are de and ctm in the first place
    n_ctm_de = length(ctm_de) 
    # count number of such genes that end up in this module
    ctm_de_module = intersect(ctm_de, ma[ma[paste0(tissue, '_module')] == module,]$gene_name) %>% unique(.)
    n_ctm_de_module = length(ctm_de_module)
    # calculate enrichment ratio
    ratio_ctm_de = (n_ctm_de_module / nmodule) / (n_ctm_de / get(paste0(tissuecode, '.nGenes')))
    # test for enrichment
    pval_ctm_de_enrich = phyper(n_ctm_de_module, n_ctm_de, get(paste0(tissuecode, '.nGenes'))-n_ctm_de, nmodule, lower.tail = F)

    if (coolmod == T) {
    
      # investigate connection between +/- ME correlation and ctm genes
      mm = get(paste0(tissuecode, '.geneModuleMembership'))[paste0('MM_', toupper(tissuecode), '_', module)]
      mm['gene_id'] = rownames(mm)
      # mm = merge(mm, ma[,c('gene_id', 'gene_name')], by = 'gene_id', all.x = T)
      mm = mm[mm$gene_id %in% modulegenes,]

      mm.pos = mm[mm[1] > 0,]$gene_id
      mm.neg = mm[mm[1] < 0,]$gene_id
      write.table(mm.pos, paste0('../../data/pos_corr_genes_', tissue, '_', module, '.txt'), quote = F, row.names = F, col.names = F)
      write.table(mm.neg, paste0('../../data/neg_corr_genes_', tissue, '_', module, '.txt'), quote = F, row.names = F, col.names = F)

      n_pos_ctm = sum(ma[ma$gene_name %in% modulectm,]$gene_id %in% mm.pos)
      n_neg_ctm = length(modulectm) - n_pos_ctm

      posneg_result = data.frame(tissue = tissue,
                                  module = module,
                                  celltype = c,
                                  n_pos_corr = length(mm.pos),
                                  n_neg_corr = length(mm.neg),
                                  n_ctm_in_module = nmodulectm,
                                  n_ctm_pos_corr = n_pos_ctm,
                                  n_ctm_neg_corr = n_neg_ctm
                                  )
      write.table(posneg_result, paste0('../../data/posneg_MEcorr_counts_', tissue, '_', module, '_', c, '.txt'), row.names = F)
    
      write.table(modulectm, paste0('../../data/ctm_overlap_', tissue, '_', module, '_', c, '.txt'), quote = F, row.names = F, col.names = F)
      write.table(modulegenes_sym[!(modulegenes_sym %in% modulectm)], paste0('../../data/NON_ctm_overlap_', tissue, '_', module, '_', c, '.txt'), quote = F, row.names = F, col.names = F)
      write.table(ctm_de_module, paste0('../../data/ctm_de_overlap_', tissue, '_', module, '_', c, '.txt'), quote = F, row.names = F, col.names = F)

    }

    thismodule_enrich = rbind(thismodule_enrich, data.frame(
                        tissue = tissue,
                        module = module,
                        celltype = c,
                        n_ctm_in_module = nmodulectm,
                        n_ctm = nctm,
                        n_in_module = nmodule,
                        ctm_enrich_ratio = ratio_ctm,
                        pval_ctm_enrich = pval_ctm_enrich,
                        n_ctm_and_de = n_ctm_de,
                        n_ctm_de_inmodule = n_ctm_de_module,
                        ctm_de_enrich_ratio = ratio_ctm_de,
                        pval_ctm_de_enrich = pval_ctm_de_enrich))

  }
  return(thismodule_enrich)
}

a.ctm_enrich = data.frame()
for (m in a.modNames) {
  a.ctm_enrich = rbind(a.ctm_enrich, test_ctm_enrichment('a', m))
}
a.ctm_enrich['pval_ctm_enrich'] = as.numeric(a.ctm_enrich$pval_ctm_enrich)
write.table(a.ctm_enrich[a.ctm_enrich$pval_ctm_enrich < 0.05,] %>% arrange(pval_ctm_enrich), '../../data/adipose_ctm_enrichment.txt', row.names = F)

l.ctm_enrich = data.frame()
for (m in l.modNames) {
  l.ctm_enrich = rbind(l.ctm_enrich, test_ctm_enrichment('l', m))
}
l.ctm_enrich['pval_ctm_enrich'] = as.numeric(l.ctm_enrich$pval_ctm_enrich)
l.ctm_enrich[l.ctm_enrich$module == 'saddlebrown' & 
            l.ctm_enrich$pval_ctm_enrich < 0.05 & 
            l.ctm_enrich$n_ctm_in_module > 0,]
write.table(l.ctm_enrich[l.ctm_enrich$pval_ctm_enrich < 0.05,] %>% arrange(pval_ctm_enrich), '../../data/liver_ctm_enrichment.txt', row.names = F)







# # compute gene significance based on cross-tissue correlation
# # iterate over all pairwise cross-tissue combinations of genes/MEs
# # adipose gene expression correlated to liver MEs
# xtissueSignif_A2L = data.frame(matrix(NA, nrow = a.nGenes, ncol = 0))
# for (liverMEidx in seq(1, ncol(l.MEs))) {
#   # select a liver ME
#   thisliverME = l.MEs[liverMEidx]
#   # compute associations to this liver module
#   geneTraitSignificance = as.data.frame(cor(a.datExpr, thisliverME, use = "p"))
#   names(geneTraitSignificance) = paste("GS.", names(thisliverME), sep="")
#   GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
#   names(GSPvalue) = paste("p.GS.", names(thisliverME), sep="")
#   xtissueSignif_A2L = cbind(xtissueSignif_A2L, geneTraitSignificance, GSPvalue)
# }
# # liver gene expression correlated to adipose MEs
# xtissueSignif_L2A = data.frame(matrix(NA, nrow = l.nGenes, ncol = 0))
# for (adiMEidx in seq(1, ncol(a.MEs))) {
#   # select an adipose ME
#   thisadiME = a.MEs[adiMEidx]
#   # compute associations to this adipose module
#   geneTraitSignificance = as.data.frame(cor(l.datExpr, thisadiME, use = "p"))
#   names(geneTraitSignificance) = paste("GS.", names(thisadiME), sep="")
#   GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
#   names(GSPvalue) = paste("p.GS.", names(thisadiME), sep="")
#   xtissueSignif_L2A = cbind(xtissueSignif_L2A, geneTraitSignificance, GSPvalue)
# }

# # identify 'key genes' in correlated modules
# MEMBERTHRESH = 0.75
# get_top_quantile <- function(vec) {
#   return( vec > quantile(vec, 0.9) | vec < quantile(vec, 0.1) )
# }
# # adipose
# a.himember = a.geneModuleMembership > MEMBERTHRESH
# a.topquant = xtissueSignif_A2L %>%  
#               select(starts_with("GS")) %>% 
#               summarise_at(vars(starts_with("GS")), get_top_quantile)
# # liver
# l.himember = l.geneModuleMembership > MEMBERTHRESH
# l.topquant = xtissueSignif_L2A %>%  
#               select(starts_with("GS")) %>% 
#               summarise_at(vars(starts_with("GS")), get_top_quantile)





#(from within top loop)
  # # generate plots showing module membership vs. cross-tissue significance
  # a.keygenes = a.himember[,paste0("MM_A_", a.xtissuemodule)] & a.topquant[,paste0("GS.L_ME", l.xtissuemodule)]
  # a_to_l_scatter = ggplot(
  #   mapping = aes(x = a.geneModuleMembership[,paste0("MM_A_", a.xtissuemodule)],
  #                 y = xtissueSignif_A2L[,paste0("GS.L_ME", l.xtissuemodule)],
  #                 color = a.keygenes)) +
  #       geom_point() + 
  #       xlab(paste0("Module membership to A_", a.xtissuemodule)) + 
  #       ylab(paste0("Correlation to L_", l.xtissuemodule))
  # ggsave(paste0("../../fig/A", a.xtissuemodule, "_to_L", l.xtissuemodule, "_scatter.png"), a_to_l_scatter)

  # l.keygenes = l.himember[,paste0("MM_L_", l.xtissuemodule)] & l.topquant[,paste0("GS.A_ME", a.xtissuemodule)]
  # l_to_a_scatter = ggplot(
  #   mapping = aes(x = l.geneModuleMembership[,paste0("MM_L_", l.xtissuemodule)],
  #                 y = xtissueSignif_L2A[,paste0("GS.A_ME", a.xtissuemodule)],
  #                 color = l.keygenes)) +
  #       geom_point() +
  #       xlab(paste0("Module membership to ", l.xtissuemodule)) + 
  #       ylab(paste0("Correlation to A_", a.xtissuemodule))
  # ggsave(paste0("../../fig/L", l.xtissuemodule, "_to_A", a.xtissuemodule, "_scatter.png"), l_to_a_scatter)

  # # overlap with gene names
  # a.keygenes_df = annot[annot$gene_id %in% names(a.keygenes[a.keygenes==T]),c("gene_name", "gene_id")]
  # l.keygenes_df = annot[annot$gene_id %in% names(l.keygenes[l.keygenes==T]),c("gene_name", "gene_id")]
  # write.table(a.keygenes_df, paste0("../../data/adipose_keygenes_A_", a.xtissuemodule, "_to_L_", l.xtissuemodule, ".txt"), quote = F, row.names = F)
  # write.table(l.keygenes_df, paste0("../../data/liver_keygenes_L_", l.xtissuemodule, "_to_A", a.xtissuemodule, ".txt"), quote = F, row.names = F)

