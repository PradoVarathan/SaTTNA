ClusterCoeffGen = function(input_data){
  N_deg_data = data.frame("gene_id" = character(),"Neighbour_connectivity" = integer(), "Degree" = integer())

  for (i in 1:nrow(input_data)) {
    count_deg = 0
    comp_row = as.data.frame(input_data[i,])
    degree = sum(comp_row)
    if (degree != 0) {
      for (k in 1:nrow(input_data)) {
        if (i != k) {
          if (input_data[i,k] == 1) {
            check_row = as.data.frame(input_data[k,])
            logvec = ((comp_row == check_row) & (comp_row == 1))
            count_deg = count_deg + sum(logvec)
          }
        }
      }
    }
    N_deg_data = rbind(N_deg_data,data.frame("gene_id" = rownames(input_data)[i],"Neighbour_connectivity" = count_deg ,"Degree" = degree))
  }

  N_deg_data$ClusteringCoeff = (2*(N_deg_data$Neighbour_connectivity))/((N_deg_data$Degree)*(N_deg_data$Degree - 1))

  return(N_deg_data)
}
