David_Analysis = function(email,GeneList){
  david = DAVIDWebService$new(email="ppugale@iu.edu",
                              url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

  mart <- useMart(biomart ="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

  attributes=c('ensembl_gene_id','hgnc_symbol','entrezgene_id')

  converted_g_list<- getBM(attributes=attributes, filters="hgnc_symbol",values=GeneList,
                           mart=mart, uniqueRows=T)

  genelist = converted_g_list$ensembl_gene_id

  result = addList(david,genelist, idType = "ENSEMBL_GENE_ID",
                   listName = "genelist", listType="Gene")

  setAnnotationCategories(david,c("GOTERM_BP_ALL"))

  termCluster = getClusterReport(david,type = "Term")

  getClusterReportFile(david,type = "Term", fileName = "termClusterReport.tab")

  All_Cluster = data.frame(termCluster@cluster[[1]]$Members)

  new_plot = plot2D(termCluster,1)

  return(new_plot)
}
