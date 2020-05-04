Reactome_Analysis = function(GeneList, type_plt = "emap"){

  mart <- useMart(biomart ="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

  attributes=c('ensembl_gene_id','hgnc_symbol','entrezgene_id')

  converted_g_list<- getBM(attributes=attributes, filters="hgnc_symbol",values=GeneList,
                           mart=mart, uniqueRows=T)
  genelist = converted_g_list$entrezgene_id
  de = as.character(genelist)
  x <- enrichPathway(gene=de,pvalueCutoff=0.01, readable=T)
  res = as.data.frame(x)
  if (type_plt == "bar"){
    plt = barplot(x,showCategory = 15)
  } else if (type_plt == "dot") {
    plt = dotplot(x, showCategory=15)
  } else if (type_plt == "emap"){
    plt = emapplot(x)
  } else if (type_plt == "cnet"){
    cnetplot(x, categorySize="pvalue", foldChange=genelist)
  }

  return(plt)

}
