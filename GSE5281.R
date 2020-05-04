#Unpack the CEL files
library(R.utils)
library(affy)
library(gcrma)
library(readxl)
library(dplyr)
library(GOSemSim)
library("illuminaHumanv4.db")
library(tibble)
library(biomaRt)
library(ggplot2)
library(gganimate)
library(ReactomePA)

setwd("/Volumes/PradoHD/GEODatasets/GSE5281_RAW/")
#data = getwd()
#cels = list.files( pattern = "CEL")
#sapply(cels, gunzip)

cels = list.files(pattern = "CEL")


raw.data = ReadAffy(verbose = TRUE, filenames = cels, cdfname = "hgu133plus2") #From bioconductor

#perform RMA normalization (I would normally use GCRMA but it did not work with this chip)
data.rma.norm = rma(raw.data)

#Get the important stuff out of the data - the expression estimates for each array
rma = exprs(data.rma.norm)

#Format values to 5 decimal places
rma = format(rma, digits = 5)

# ON PIPE -----------------------------------------------------------------

# Preparing input data ----------------------------------------------------

gse5281 = as.data.frame(rma)
rm(rma,raw.data,data.rma.norm)

# Converting affyIDs to entrez IDs and averaging genes if dupl --------


gse5281$affy_hg_u133_plus_2 = rownames(gse5281)
gse5281 = na.omit(gse528)
typeof(gse5281$GSM119615.CEL)
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(mart = mart, attributes = c("affy_hg_u133_plus_2", "entrezgene_id"), filter = "affy_hg_u133_plus_2",
                     values = rownames(gse5281),uniqueRows = TRUE)

gse5281 = merge(gse5281, annotLookup, by = 'affy_hg_u133_plus_2')
gse5281 = na.omit(gse5281)
gse5281$affy_hg_u133_plus_2 = NULL
for (i in 1:length(gse5281)) {
  
  if ( colnames(gse5281)[i] != "entrezgene_id") {
    
    gse5281[,i] = as.numeric(as.character(gse5281[,i] ))
    
  }
  
}
gse5281 = as.data.frame(gse5281)
gse5281 = gse5281 %>% group_by(entrezgene_id) %>% mutate_each(funs(mean)) %>% distinct

gse5281 = data.frame(column_to_rownames(gse5281, var = "entrezgene_id"))
gse5281$entrezgene_id = NULL

#rm(mart,sample,annotLookup)






# Calculating variances over samples for each gene ------------------------

gse5281_variance = GeneSampleVariance(gse5281)

# Vairance plots threshold visualization ----------------------------------

variance_plot = Gene_Variance_Plot(gse5281_variance, start = 1000)
#variance_plot
#anim_save("Variance_gse5281.gif")
#chose 0.9
rm(variance_plot)

# Prepare data with variance and sample information -----------------------
#setwd("/Volumes/PradoHD/GEODatasets/GSE5281_RAW/")
subject_data = read_xls("GSE5281_sample_characteristics.xls")
subject_data$Primary_Key = "nothinfornow"
for (i in 1:nrow(subject_data)) {
  
  subject_data$Primary_Key[i] = toString(c(as.character(subject_data$'Organ Region:'[i]),as.character(subject_data$'Disease State:'[i])),sep=",")
  subject_data$'GEO Accession:'[i] = paste(as.character(subject_data$'GEO Accession:'[i]), ".CEL", sep = "")
}

subject_data = data.frame("sample_id" = subject_data$'GEO Accession:' , "condition" = subject_data$Primary_Key)
colnames(subject_data) = c("sample_id","condition")
gse5281_filtered_ND = Var_Sam_Preparation(gse5281,gse5281_variance, subject_data = subject_data, for_subject = as.character(subject_data$condition[14]), threshold = 0.9)
gse5281_filtered_AD = Var_Sam_Preparation(gse5281,gse5281_variance, subject_data = subject_data, for_subject = as.character(subject_data$condition[86]), threshold = 0.9)

# Creating network based on pearson correlation ---------------------------

gse5281_network_AD = GenNet_Cor(gse5281_filtered_AD)
gse5281_network_ND = GenNet_Cor(gse5281_filtered_ND)

# Getting the semantic score matrix from GOSemSim -------------------------

semantic_score_MF_AD = GetSem_Score(gse5281_network_AD, type = "MF")
semantic_score_MF_ND = GetSem_Score(gse5281_network_ND, type = "MF")

# Preparing files - CLustering Coeffecient  & Neighbour lists -------------

gse5281_cc_AD = ClusterCoeffGen(gse5281_network_AD)
gse5281_cc_ND = ClusterCoeffGen(gse5281_network_ND)

gse5281_neighbours_AD = Extract_neighbours(gse5281_network_AD)
gse5281_neighbours_ND = Extract_neighbours(gse5281_network_ND)

# Creating modules based on the network -----------------------------------

gse5281_mods_AD = Extract_Mods(network_file = gse5281_network_AD,sem_score = semantic_score_MF_AD,clusteringcoeff = gse5281_cc_AD,neighbours_file = gse5281_neighbours_AD,
                              SST = 0.5, CCT = 0.5)
gse5281_mods_AD = Extract_Mods(network_file = gse5281_network_AD,sem_score = semantic_score_MF_AD,clusteringcoeff = gse5281_cc_AD,neighbours_file = gse5281_neighbours_AD,
                               SST = 0.5, CCT = 0.5)


