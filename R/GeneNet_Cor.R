GenNet_Cor = function(data,method=c("pearson"), threshold = 0.5, cor_req = T){
  if (cor_req == T) {
    data = cor(as.data.frame(data),method = method)
  }
  #creatin network
  output_network = data
  output_network[output_network >= threshold] = 1
  output_network[output_network < threshold] = 0

  return(output_network)
}
