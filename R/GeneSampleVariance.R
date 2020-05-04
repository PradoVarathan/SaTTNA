GeneSampleVariance = function(input_data){

  tot_genes = nrow(input_data)
  tot_subjs = length(input_data)
  variance_dataframe = data.frame("GeneID" = character(), "Variance" = numeric(),stringsAsFactors = F)

  for (g in 1:tot_genes) {
    variance = var(t(input_data[g,]))[1,1]
    geneid = rownames(input_data[g,])
    variance_dataframe[g,1] = geneid
    variance_dataframe[g,2] = variance
  }
  input_data
  # find a faster way
  return(variance_dataframe)

}
