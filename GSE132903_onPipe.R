setwd("/Volumes/PradoHD/GEODatasets/")

library(readxl)
library(dplyr)
library(GOSemSim)
library("illuminaHumanv4.db")
library(tibble)
library(biomaRt)
library(ggplot2)
library(gganimate)

# Preparing input data ----------------------------------------------------

setwd("/Volumes/PradoHD/GEODatasets/")
#list.files()
#GSE118553 = read.delim("GSE118553_non-normalized_data.txt",header = T, sep = "\t")
GSE132903 = read_xlsx("GSE132903_Matrix_Normalized.xlsx")  #middle temporal gyrus
GSE132903$ID_REF =  NULL
colnames(GSE132903)[1] = "Illumina_Probe_Name"

# Converting illuminaIDs to entrez IDs and averaging genes if dupl --------

gse132903 = Illumina_Enterez(GSE132903)
rm(GSE132903)

# Calculating variances over samples for each gene ------------------------

gse132903_variance = GeneSampleVariance(gse132903)

# Vairance plots threshold visualization ----------------------------------

variance_plot = Gene_Variance_Plot(gse132903_variance, start = 10000)
variance_plot
#anim_save("Variance_GSE132903.gif")

# Prepare data with variance and sample information -----------------------

subject_data = data.frame("sample_id" = as.data.frame(colnames(gse132903)), "condition" = as.data.frame(as.character(getColNames(gse132903,char = "[_]", pos = 1))))
colnames(subject_data) = c("sample_id","condition")
gse132903_filtered_AD = Var_Sam_Preparation(gse132903,gse132903_variance, subject_data = subject_data, for_subject = "AD", threshold = 0.5)
gse132903_filtered_ND = Var_Sam_Preparation(gse132903,gse132903_variance, subject_data = subject_data, for_subject = "ND", threshold = 0.5)

# Creating network based on pearson correlation ---------------------------

gse132903_network_AD = GenNet_Cor(gse132903_filtered_AD)
gse132903_network_ND = GenNet_Cor(gse132903_filtered_ND)

# Getting the semantic score matrix from GOSemSim -------------------------

semantic_score_MF_AD = GetSem_Score(gse132903_network_AD, type = "MF")
semantic_score_MF_ND = GetSem_Score(gse132903_network_ND, type = "MF")

# Preparing files - CLustering Coeffecient  & Neighbour lists -------------

gse132903_cc_AD = ClusterCoeffGen(gse132903_network_AD)
gse132903_cc_ND = ClusterCoeffGen(gse132903_network_ND)

gse132903_neighbours_AD = Extract_neighbours(gse132903_network_AD)
gse132903_neighbours_ND = Extract_neighbours(gse132903_network_ND)

# Creating modules based on the network -----------------------------------

# Modules_AD = Extract_Mods(gse132903_network_AD,semantic_score_MF_AD,gse132903_neighbours_AD,gse132903_cc_AD)
# Modules_ND = Extract_Mods(gse132903_network_ND,semantic_score_MF_ND,gse132903_neighbours_ND,gse132903_cc_ND)

# Showing the modules and converting them to GeneID -----------------------



# Getting Reactome Networks for top module --------------------------------



# Analysis with DAVID functional tool for Top module ----------------------



