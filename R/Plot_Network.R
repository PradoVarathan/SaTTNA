Plot_Network = function(Modules,network_file){
  Modules = Modules
  network_Cnt = network_file

  genes_mod_1 = as.character(Modules[1,2])
  genes_mod_1 = trimws(unlist(strsplit(genes_mod_1,",")))
  genes_mod_1 = unique(genes_mod_1)
  gse132903_mod_1 = as.data.frame(network_Cnt[genes_mod_1,genes_mod_1])

  genelist = colnames(gse132903_mod_1)
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  attributes = c('hgnc_symbol','entrezgene_id')
  converted_g_list <- getBM(attributes = attributes, filters = "entrezgene_id",values = genelist,
                            mart = mart)

  GeneList_Cnt = converted_g_list$hgnc_symbol

  converted_g_list = converted_g_list %>% group_by(entrezgene_id) %>% mutate(hgnc_symbol= paste0(hgnc_symbol, collapse = ","))
  converted_g_list = converted_g_list[!duplicated(converted_g_list$entrezgene_id),]
  gse132903_mod_1$entrezgene_id = rownames(gse132903_mod_1)
  gse132903_mod_1 = merge(gse132903_mod_1,converted_g_list, by = "entrezgene_id")

  rownames(gse132903_mod_1) = gse132903_mod_1$hgnc_symbol
  gse132903_mod_1$entrezgene_id = gse132903_mod_1$hgnc_symbol = NULL

  gse132903_mod_1 = as.data.frame(t(gse132903_mod_1))

  gse132903_mod_1$entrezgene_id = rownames(gse132903_mod_1)
  gse132903_mod_1 = merge(gse132903_mod_1,converted_g_list, by = "entrezgene_id")
  rownames(gse132903_mod_1) = gse132903_mod_1$hgnc_symbol
  gse132903_mod_1$entrezgene_id = gse132903_mod_1$hgnc_symbol = NULL




  # Creating the network plots ----------------------------------------------

  mod1 <- graph_from_adjacency_matrix(as.matrix(gse132903_mod_1),)
  # plot(mod1,vertex.size=1,edge.arrow.size = 0.5,edge.arrow.width = 0.1,layout = layout_with_dh, main ="Module 1",edge.width=0.3)
  mod1 <- simplify(mod1, remove.multiple = F, remove.loops = T)
  return(mod1)

}
