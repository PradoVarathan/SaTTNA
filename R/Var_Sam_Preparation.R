
Var_Sam_Preparation = function(input_data, input_var_data,subject_data, threshold = 0.5, for_subject = character()){

  input_data$GeneID = rownames(input_data)
  input_data = merge(input_data,input_var_data, by = "GeneID")
  input_data = filter(input_data, Variance > threshold)
  rownames(input_data) = input_data$GeneID

  #input_data = data.frame(column_to_rownames(input_data, var = "GeneID"))
  input_data$GeneID = NULL
  input_data = as.data.frame(t(input_data))
  input_data$sample_id = rownames(input_data)
  input_data = merge(input_data,subject_data,by = "sample_id")
  input_data$condition = as.character(input_data$condition)
  input_data_subj = filter(input_data,condition == for_subject)
  rownames(input_data_subj) = input_data_subj$sample_id
  input_data_subj$sample_id = input_data_subj$condition = NULL

  return(input_data_subj)
}

