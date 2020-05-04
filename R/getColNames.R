getColNames = function(data , char = "X", pos = 2, start = 1){
  k = start
  out = c()
  for (i in k:length(data)) {
    out = c(out,unlist(strsplit(as.character(colnames(data)[i]),char))[[pos]])
    k = k + 1
  }

  return(out)
}
