## Pseudoinverse calculations on weighted matrices on heavy networks






rm(list = ls())
if (require("rstudioapi")) {
  script_path <- rstudioapi::getSourceEditorContext()$path
}

library(tidyr)
library(igraph)
library(vegan)
library(dplyr)
library(purrr)
library(corpcor)
library(Matrix)
source("null_model_generation_functions.R")

script_directory <- dirname(script_path)
setwd(script_directory)
options(max.print=10000,width=100)

decompose_matrix <- function(original_matrix) {
  # Create a binary matrix by setting non-zero values to 1
  column_names <- colnames(original_matrix)
  row_names <- rownames(original_matrix)
  new_column_names <- c(row_names,column_names)
  #print(new_column_names)
  binary_matrix <- matrix(0, nrow = (nrow(original_matrix)*ncol(original_matrix)), ncol = ncol(original_matrix)+nrow(original_matrix),dimnames = list(NULL,new_column_names))
  #print(binary_matrix)
  log_vector <- c()
  aux_count <- 0
  for (i in 1:nrow(original_matrix)) {
    for (j in 1:ncol(original_matrix)) {
      aux_count <- aux_count + 1
      if(original_matrix[i, j]==0)
      {
        log_vector <- c(log_vector,0)
      }
      else
      {
        log_vector <- c(log_vector, log(original_matrix[i, j]))
        col_index_1 <- match(row_names[i],new_column_names)
        col_index_2 <- match(column_names[j],new_column_names)
        binary_matrix[aux_count,col_index_1] <- 1
        binary_matrix[aux_count,col_index_2] <- 1
      }
      
      
    }
  }
  #print(binary_matrix[binary_matrix!=0])
  #sparse_matrix <- sparseMatrix(i = which(binary_matrix != 0), 
                                #j = which(binary_matrix != 0),
                                #x = binary_matrix[binary_matrix != 0])
  #print(sparse_matrix)
  #print(sparse_matrix[sparse_matrix!=0])
  return(list(binary_matrix = binary_matrix, log_vector = log_vector))
}

calculate_final_interaction_matrix <- function(interaction_matrix,effective_abundance_vector)
{
  final_matrix <- matrix(0, nrow = nrow(interaction_matrix), ncol = ncol(interaction_matrix),dimnames = list(rownames(interaction_matrix),colnames(interaction_matrix)))
  for (i in 1:nrow(interaction_matrix)) {
    for (j in 1:ncol(interaction_matrix)) {
      final_matrix[i,j] <- interaction_matrix[i,j]/(effective_abundance_vector[i]*effective_abundance_vector[j])
    }
  }
  return(final_matrix)
}



datadir <- "data/"
csv_directory <- "./data"

# Get a list of all .csv files in the directory
csv_files <- list.files(csv_directory, pattern = "\\.csv$", full.names = FALSE)
csv_files

# Remove the .csv extension and store the file names in a vector

file_names <- sub("\\.csv$", "", csv_files) 

# Print the vector of file names
print(file_names)
## If you want to specify a model, substitute this following line 
## by a vector of models locations in file system as shown below
lnetw <- file_names

for (netw in (lnetw)){
  
  print(netw)
  result_analysis <- analyze_network(paste(netw,".csv",sep=""), directory = datadir, guild_a = "Plant", guild_b = "Pollinator", only_NODF = TRUE)
  weighted_network <- sum(result_analysis$matrix > 1)>0
  if(weighted_network)
  {
    result <- decompose_matrix(result_analysis$matrix)
    roger_penrose_inverse <- pseudoinverse(result$binary_matrix)
    #print(class(roger_penrose_inverse))
    #print(result$log_vector)
    effective_abundance <- roger_penrose_inverse  %*% result$log_vector
    print(effective_abundance)
    #print(effective_abundance)
    #print(result_analysis$matrix)
    #print("This is matrix calculated as said in the paper")
    #print(effective_abundance)
    
    final_mass_action_matrix <- calculate_final_interaction_matrix(result_analysis$matrix,effective_abundance)
    write.csv(final_mass_action_matrix, file = paste0(paste0("interaction_abundance/",netw),".csv"))
    #print(final_mass_action_matrix)
    #print("This is final mass-action matrix accounting for effective abundance")
    #print(final_mass_action_matrix)
  }
}