Illumina_Enterez = function(input_data){

  ProbeIds = as.character(input_data$Illumina_Probe_Name)
  Gene_Ids = data.frame(Gene = unlist(mget(x = ProbeIds,envir = illuminaHumanv4SYMBOL)), "Illumina_Probe_Name" = ProbeIds)
  input_data = merge(input_data, Gene_Ids, by = "Illumina_Probe_Name")
  input_data$Illumina_Probe_Name = NULL
  input_data = na.omit(input_data)

  genelist = input_data$Gene
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  attributes = c('hgnc_symbol','entrezgene_id')
  converted_g_list <- getBM(attributes = attributes, filters = "hgnc_symbol",values = genelist,
                            mart = mart, uniqueRows = T)
  colnames(converted_g_list)[1] = "Gene"

  input_data = merge(input_data, converted_g_list, by = "Gene")
  input_data$Gene = NULL
  input_data = na.omit(input_data)
  input_data = input_data %>% group_by(entrezgene_id) %>% summarise_each(funs(mean))
  input_data = data.frame(column_to_rownames(input_data, var = "entrezgene_id"))
  input_data$entrezgene_id = NULL
  return(input_data)


}
