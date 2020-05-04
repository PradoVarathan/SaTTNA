GetSem_Score = function(input_data, type = "MF") {
  gene_list = colnames(input_data)
  semData = godata('org.Hs.eg.db',ont = type, computeIC = FALSE)
  Sem_Score = mgeneSim(gene_list,semData = semData,measure = "Wang")
  return(Sem_Score)
}
