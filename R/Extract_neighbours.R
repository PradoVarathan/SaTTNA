Extract_neighbours = function(network_file) {

  neighbours_file = data.frame("gene_id" = character(),"Neighbours" = list())
  for (i in 1:nrow(network_file)) {
    g_s = rownames(network_file)[i]
    n_s = colnames(network_file)[as.logical(network_file[i,])]
    n_s = n_s[ n_s != g_s]
    neighbours_file = rbind(neighbours_file,data.frame("gene_id" = g_s,
                                                       "Neighbours" = toString(n_s,sep = ",")))

  }
  return(neighbours_file)
}
