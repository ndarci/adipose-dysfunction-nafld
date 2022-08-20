module = commandArgs(trailingOnly = T)[1]

# module = 'A_lightyellow'

# import panther results
panth = read.table(paste0('../../data/pantherGeneList_', module, '.txt'), sep = '\t', fill = T)
names(panth) = c('gene_id', 'mapped_ids', 'name_symbol_ortholog', 'panther_family', 'panther_prot_class', 'species')

# remember which module these genes came from
panth['module'] = module

# extract gene symbol
splname = strsplit(panth$name_symbol_ortholog, ';', fixed = T)
panth['gene_symbol'] = sapply(splname, '[[', 2)

# output TFs
write.table(panth[grepl('transcription factor', panth$panther_prot_class),c('module', 'gene_symbol', 'panther_prot_class')], paste0("../../data/TFs_", module, ".txt"))
