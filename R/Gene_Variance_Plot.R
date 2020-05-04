Gene_Variance_Plot = function(input_data,start = 5000){
  var_0.4 = filter(input_data, Variance > 0.4)
  var_0.6 = filter(input_data, Variance > 0.6)
  var_0.8 = filter(input_data, Variance > 0.8)
  var_1.0 = filter(input_data, Variance > 1)
  var_1.2 = filter(input_data, Variance > 1.2)
  var_1.4 = filter(input_data, Variance > 1.4)
  var_1.6 = filter(input_data, Variance > 1.6)
  var_1.8 = filter(input_data, Variance > 1.8)
  var_2.0 = filter(input_data, Variance > 2)


  a <- data.frame(Variances = c("0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0"),
                  No_of_genes = c(start,start,start,start,start,start,start,start,start),
                  frame = rep('a',9))
  b <- data.frame(Variances = c("0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0"),
                  No_of_genes = c(nrow(var_0.4),nrow(var_0.6),nrow(var_0.8),nrow(var_1.0),nrow(var_1.2),nrow(var_1.4),nrow(var_1.6),nrow(var_1.8),nrow(var_2.0)),
                  frame = rep('b',9))
  data <- rbind(a,b)
  ggplot(a, aes(x = Variances, y = No_of_genes, fill = Variances)) +
    geom_bar(stat = 'identity')
  plot_var = ggplot(data, aes(x = Variances, y = No_of_genes, fill = Variances)) +
    geom_bar(stat = 'identity') +
    theme_bw() + transition_states(frame,transition_length = 2,state_length = 1) + ease_aes('sine-in-out')

  return(plot_var)

}
