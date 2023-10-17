
###### Considerations ######
#As there are many ways to handle the calculations of the normalized data, they are all included below
#The latest results corresponds to the execution of the function in the loop for that are not commented
#In case there is a need to calculate them in other way, the methods are still provided

##### Reset memory and charge libraries ####
rm(list = ls())
if (require("rstudioapi")) {
  script_path <- rstudioapi::getSourceEditorContext()$path
}
#options(max.print=2000)
script_directory <- dirname(script_path)
setwd(script_directory)
library(tidyr)
library(igraph)
library(vegan)
library(dplyr)
library(purrr)
library(ggrepel)


##### Functions #####

source("null_model_generation_functions.R")


###################
#Global lists that keep data to be processed in different ways
#- This is not scalable, but it is sufficient to satisfy the number of networks available
###################
data_lists <- list()
data_lists_all_graphs <- list()
data_lists_all_graphs_laplacian <- list()
data_lists_all_graphs_laplacian_normalized <- list()
data_lists_laplacian <- list()
data_lists_laplacian_normalized <- list()
# Function to add data to the named list 
add_data <- function(name, data) {
  if (!(name %in% names(data_lists))) {
    data_lists[[name]] <<- c()
  }
  data_lists[[name]] <<- c(data_lists[[name]], data)
}
add_data_laplacian <- function(name, data) {
  if (!(name %in% names(data_lists_laplacian))) {
    data_lists_laplacian[[name]] <<- c()
  }
  data_lists_laplacian[[name]] <<- c(data_lists_laplacian[[name]], data)
}
add_data_laplacian_normalized <- function(name, data) {
  if (!(name %in% names(data_lists_laplacian_normalized))) {
    data_lists_laplacian_normalized[[name]] <<- c()
  }
  data_lists_laplacian_normalized[[name]] <<- c(data_lists_laplacian_normalized[[name]], data)
}

add_data_all_graphs <- function(name, list_to_be_added,data)
{
  if (!(name %in% names(list_to_be_added))) {
    list_to_be_added[[name]] <- c()
  }
  list_to_be_added[[name]] <- c(list_to_be_added[[name]], data)
  return(list_to_be_added)
}

list_of_normalization_criterias <- list()
add_list <- function(outer_name, inner_name, values) {
  if (!(outer_name %in% names(list_of_normalization_criterias))) {
    list_of_normalization_criterias[[outer_name]] <<- list()
  }
  
  inner_list <- list(name = inner_name, values = values)
  list_of_normalization_criterias[[outer_name]] <<- c(list_of_normalization_criterias[[outer_name]], inner_list)
}

## MIN-MAX NORMALIZATION CRITERIA FUNCTIONS ####

min_max_normalization <- function(plot_filename="",data_lists,network)
{
  final_distances_list <- list()
  print("######## MIN-MAX NORMALIZATION DATA #############")
  for (name in names(data_lists)) {
    cat("Data for name:", name, "\n")
    cat("Spectral distances calculated:",data_lists[[name]],"\n\n")
    min_distance_found <-  min(data_lists[[name]])
    max_distance_found <- max(data_lists[[name]])
    final_distances <- c()
    for(distance in data_lists[[name]]){
      new_distance <- (distance -min_distance_found)/(max_distance_found-min_distance_found)
      final_distances <- c(final_distances,new_distance)
      
    }
    final_distances_list[[name]] <- final_distances
    cat("Spectral distances normalized:",final_distances,"\n\n")
  }
  plot_box_diagrams_normalized(paste0(network,"-MinMax-Criteria"),final_distances_list,plot_filename)
}
min_max_normalization_model_independant <- function(plot_filename="",data_lists,network)
{
  min_value <- Inf
  max_value <- 0
  final_distances_list <- list()
  print("######## MIN-MAX NORMALIZATION DATA #############")
  for (name in names(data_lists)) {
    vector_min <- min(data_lists[[name]])# Find the minimum value in the current vector
    vector_max <- max(data_lists[[name]])
    if (vector_min < min_value) {
      min_value <- vector_min  # Update the minimum value if a smaller value is found
    }
    if(vector_max > max_value)
    {
      max_value <- vector_max
    }
  }
  for (name in names(data_lists)) {
    cat("Data for name:", name, "\n")
    cat("Spectral distances calculated:",data_lists[[name]],"\n\n")
    final_distances <- c()
    for(distance in data_lists[[name]]){
      new_distance <- (distance -min_value)/(max_value-min_value)
      final_distances <- c(final_distances,new_distance)
      
    }
    final_distances_list[[name]] <- final_distances
    cat("Spectral distances normalized:",final_distances,"\n\n")
  }
  plot_box_diagrams_normalized(paste0(network,"-MinMax-ByNetwork-Criteria"),final_distances_list,plot_filename)
}

min_max_normalization_all_distances <- function(data_list_all_graphs_x, max_value, min_value ,directory)
{
  for(model in names(data_list_all_graphs_x))
  {
    actual_folder_path <- file.path(folder_path_list[[model]],directory)
    final_distances_list <- list()
    print("######## MIN-MAX ALL DISTANCES DATA #############")
    for(name in names(data_list_all_graphs_x[[model]]))
    {
      final_distances <- c()
      for(distance in data_list_all_graphs_x[[model]][[name]]){
        new_distance <- (distance -min_value)/(max_value-min_value)
        final_distances <- c(final_distances,new_distance)
        
      }
      final_distances_list[[name]] <- final_distances
    }
    plot_box_diagrams_normalized(paste0(model,"-MinMax-AllDistances-Criteria"),final_distances_list,actual_folder_path)
  }
}
#### 

###### Bounds Normalization Criteria Functions  ######
## These methods have been discarded by results, however they are still provided here

spectral_distances_bounds <- list()
spectral_distances_bounds_laplacian <- list()
spectral_distances_bounds_laplacian_normalized <- list()

calculate_and_save_spdbound <- function(name, bounds,spectral_distance,length_of_sum) {
  if (!(name %in% names(spectral_distances_bounds))) {
    spectral_distances_bounds[[name]] <<- c()
  }
  spectral_distance_normalized <- spectral_distance/sqrt(sum(rep((bounds$max-bounds$min)^2,length_of_sum)))
  spectral_distances_bounds[[name]] <<- c(spectral_distances_bounds[[name]], spectral_distance_normalized)
  
  
}
calculate_and_save_spdbound_laplacian <- function(name, bounds,spectral_distance,length_of_sum) {
  if (!(name %in% names(spectral_distances_bounds_laplacian))) {
    spectral_distances_bounds_laplacian[[name]] <<- c()
  }
  spectral_distance_normalized <- spectral_distance/sqrt(sum(rep((bounds$max-bounds$min)^2,length_of_sum)))
  spectral_distances_bounds_laplacian[[name]] <<- c(spectral_distances_bounds_laplacian[[name]], spectral_distance_normalized)
  
  
}
calculate_and_save_spdbound_laplacian_normalized <- function(name, bounds,spectral_distance,length_of_sum) {
  if (!(name %in% names(spectral_distances_bounds_laplacian_normalized))) {
    spectral_distances_bounds_laplacian_normalized[[name]] <<- c()
  }
  spectral_distance_normalized <- spectral_distance/sqrt(sum(rep((bounds$max-bounds$min)^2,length_of_sum)))
  spectral_distances_bounds_laplacian_normalized[[name]] <<- c(spectral_distances_bounds_laplacian_normalized[[name]], spectral_distance_normalized)
  
  
}

calculate_bounds_adjacency <- function(model_adjacency_matrix,null_adjacency_matrix)
{
  ones_vector_model <- rep(1, ncol(model_adjacency_matrix))
  ones_vector_null <- rep(1, ncol(null_adjacency_matrix))
  max_result_model <- max(model_adjacency_matrix %*% ones_vector_model)
  min_result_null <- min(null_adjacency_matrix %*% ones_vector_null)
  return(list(max=max_result_model,min=min_result_null))
}

bounds_normalization <- function(plot_filename="",network) ## show function
{
  print("######## BOUNDS NORMALIZATION DATA #############")
  for (name in names(data_lists)) {
    cat("Data for name:", name, "\n")
    cat("Spectral distances calculated:",data_lists[[name]],"\n\n")
    cat("Spectral distances normalized:",spectral_distances_bounds[[name]],"\n\n")
    
  }
  plot_box_diagrams_normalized(paste0(network,"-Bounds-Criteria"),spectral_distances_bounds,plot_filename)
}
#####

###### Bounds Normalization Criteria - Laplacian ######
### Adjacency, Laplacian and Laplacian normalized bounds are calculated in different but efficient ways according to the results of 
### "Lectures of Network Systems" by Francesco Bullo
actual_network_max_eigenvalue_laplacian <- 0
actual_network_number_eigenvalues <- 1
actual_network_max_eigenvalue_laplacian_normalized <- 0
actual_network_number_eigenvalues_normalized <- 1

save_parameters_laplacian <- function(length,value)
{
  actual_network_max_eigenvalue_laplacian <<- value
  actual_network_number_eigenvalues <<- length
}
save_parameters_laplacian_normalized <- function(length,value)
{
  actual_network_max_eigenvalue_laplacian_normalized <<- value
  actual_network_number_eigenvalues_normalized <<- length
}
bounds_normalization_laplacian <- function(plot_filename="",have_weights=TRUE,network)
{
  final_distances_list <- list()
  print("######## BOUNDS NORMALIZATION DATA LAPLACIAN #############")
  if(have_weights)
  {
    for (name in names(data_lists_laplacian)) {
      cat("Data for name:", name, "\n")
      cat("Spectral distances calculated:",data_lists_laplacian[[name]],"\n\n")
      final_distances <- c()
      
      for(distance in data_lists_laplacian[[name]]){
        new_distance <- distance/sqrt(actual_network_number_eigenvalues*(actual_network_max_eigenvalue_laplacian)^2)
        final_distances <- c(final_distances,new_distance)
      }
      final_distances_list[[name]] <- final_distances
      cat("Spectral distances normalized:",final_distances,"\n\n")
    }
    plot_box_diagrams_normalized(paste0(network,"-Bounds-Criteria"),final_distances_list,plot_filename)
  }
  else
  {
    for (name in names(data_lists_laplacian)) {
      cat("Data for name:", name, "\n")
      cat("Spectral distances calculated:",data_lists_laplacian[[name]],"\n\n")
      cat("Spectral distances normalized:",spectral_distances_bounds_laplacian[[name]],"\n\n")
      
    }
    plot_box_diagrams_normalized(paste0(network,"-Bounds-Criteria"),spectral_distances_bounds_laplacian,plot_filename)
    }
  
  
}

####

###### Bounds Normalization Criteria - Laplacian Normalized ######
bounds_normalization_laplacianNormalized <- function(plot_filename="",requirements=FALSE,have_weights=FALSE)
{
  if(have_weights)
  {
    if(requirements)
    {
      final_distances_list <- list()
      print("######## BOUNDS NORMALIZATION DATA LAPLACIAN NORMALIZED #############")
      for (name in names(data_lists_laplacian_normalized)) {
        cat("Data for name:", name, "\n")
        cat("Spectral distances calculated:",data_lists_laplacian_normalized[[name]],"\n\n")
        final_distances <- c()
        
        for(distance in data_lists_laplacian_normalized[[name]]){
          new_distance <- distance/sqrt(actual_network_number_eigenvalues_normalized*4)
          final_distances <- c(final_distances,new_distance)
        }
        final_distances_list[[name]] <- final_distances
        cat("Spectral distances normalized:",final_distances,"\n\n")
        
      }
      plot_box_diagrams_normalized("Bounds-LaplacianNormalized-Criteria",final_distances_list,plot_filename)
    }
    else
    {
      final_distances_list <- list()
      print("######## BOUNDS NORMALIZATION DATA LAPLACIAN NORMALIZED #############")
      for (name in names(data_lists_laplacian_normalized)) {
        cat("Data for name:", name, "\n")
        cat("Spectral distances calculated:",data_lists_laplacian_normalized[[name]],"\n\n")
        final_distances <- c()
        
        for(distance in data_lists_laplacian_normalized[[name]]){
          new_distance <- distance/sqrt(actual_network_number_eigenvalues_normalized*(actual_network_max_eigenvalue_laplacian_normalized)^2)
          final_distances <- c(final_distances,new_distance)
        }
        final_distances_list[[name]] <- final_distances
        cat("Spectral distances normalized:",final_distances,"\n\n")
        
      }  
    }
  }
  else
  {
    for (name in names(data_lists_laplacian_normalized)) {
      cat("Data for name:", name, "\n")
      cat("Spectral distances calculated:",data_lists_laplacian_normalized[[name]],"\n\n")
      cat("Spectral distances normalized:",spectral_distances_bounds_laplacian_normalized[[name]],"\n\n")
      
    }
    plot_box_diagrams_normalized("Bounds-Criteria",spectral_distances_bounds_laplacian_normalized,plot_filename)
  }
  
  
}

####



### IQR Normalization Functions ### 
### This part has been discarded by the latest results. However, methods are still implemented here for the record.
iqr_normalization <- function(plot_filename="",data_lists,network)
{
  final_distances_list <- list()
  print("######## IQR NORMALIZATION DATA #############")
  for (name in names(data_lists)) {
    cat("Data for name:", name, "\n")
    cat("Spectral distances calculated:",data_lists[[name]],"\n\n")
    first_quartile <- as.numeric(paste0(quantile(data_lists[[name]], probs = 0.25)))
    iqr_value <- IQR(data_lists[[name]])
    final_distances <- c()
    
    for(distance in data_lists[[name]]){
      new_distance <- (distance -first_quartile)/(iqr_value)
      final_distances <- c(final_distances,new_distance)
    }
    final_distances_list[[name]] <- final_distances
    cat("Spectral distances normalized:",final_distances,"\n\n")
  }
  plot_box_diagrams_normalized(paste0(network,"-IQR-Criteria"),final_distances_list,plot_filename)
}

iqr_normalization_all_distances_bymodeltype <- function(data_lists_all_graphs_x, directory)
{
  temporal_distances <- list()
  print("######## IQR NORMALIZATION ALL DISTANCES DATA #############")
  for(model in names(data_lists_all_graphs_x))
  {
    for(name in names(data_lists_all_graphs_x[[model]]))
    {
      for(distance in data_lists_all_graphs_x[[model]][[name]])
      {
        temporal_distances[[name]] <- c(temporal_distances[[name]],distance)
      }
      
    }
  }
  
  for(model in names(data_lists_all_graphs_x))
  {
    actual_folder_path <- file.path(folder_path_list[[model]],directory)
    final_distances_list <- list()
    for(name in names(data_lists_all_graphs_x[[model]]))
    {
      final_distances <- c()
      first_quartile <- as.numeric(paste0(quantile(temporal_distances[[name]], probs = 0.25)))
      iqr_value <- IQR(temporal_distances[[name]])
      for(distance in data_lists_all_graphs_x[[model]][[name]]){
        new_distance <- (distance -first_quartile)/(iqr_value)
        final_distances <- c(final_distances,new_distance)
      }
      final_distances_list[[name]] <- final_distances
    }
    plot_box_diagrams_normalized(paste0(model,"-IQR-AllDistances-ByModelType-Criteria"),final_distances_list,actual_folder_path)
  }
  
}
iqr_normalization_all_distances_bynetwork <- function(plot_filename="",data_lists,network)
{
  temporal_distances <- c()
  print("######## IQR NORMALIZATION ALL DISTANCES DATA #############")
  for(name in names(data_lists))
  {
    for(distance in data_lists[[name]])
    {
      temporal_distances<- c(temporal_distances,distance)
    }
    
  }
  first_quartile <- as.numeric(paste0(quantile(temporal_distances, probs = 0.25)))
  iqr_value <- IQR(temporal_distances)
  final_distances_list <- list()
  for(name in names(data_lists))
  {
    final_distances <- c()
    
    for(distance in data_lists[[name]]){
      new_distance <- (distance -first_quartile)/(iqr_value)
      final_distances <- c(final_distances,new_distance)
    }
    final_distances_list[[name]] <- final_distances
  }
  plot_box_diagrams_normalized(paste0(network,"-IQR-ByNetwork-Criteria"),final_distances_list,plot_filename)
  
}
iqr_normalization_all_distances <- function(data_lists_all_graphs_x, directory)
{
  temporal_distances <- c()
  print("######## IQR NORMALIZATION ALL DISTANCES DATA #############")
  for(model in names(data_lists_all_graphs_x))
  {
    for(name in names(data_lists_all_graphs_x[[model]]))
    {
      for(distance in data_lists_all_graphs_x[[model]][[name]])
      {
        temporal_distances<- c(temporal_distances,distance)
      }
      
    }
  }
  first_quartile <- as.numeric(paste0(quantile(temporal_distances, probs = 0.25)))
  iqr_value <- IQR(temporal_distances)
  for(model in names(data_lists_all_graphs_x))
  {
    actual_folder_path <- file.path(folder_path_list[[model]],directory)
    final_distances_list <- list()
    for(name in names(data_lists_all_graphs_x[[model]]))
    {
      final_distances <- c()
      
      for(distance in data_lists_all_graphs_x[[model]][[name]]){
        new_distance <- (distance -first_quartile)/(iqr_value)
        final_distances <- c(final_distances,new_distance)
      }
      final_distances_list[[name]] <- final_distances
    }
    plot_box_diagrams_normalized(paste0(model,"-IQR-AllDistances-Criteria"),final_distances_list,actual_folder_path)
  }
  
}
####

##### INDICES Normalization #####

indices_normalization <- function(..., plot_filename="",data_lists)
{
  final_distances_list <- list()
  print("######## INDICES NORMALIZATION DATA #############")
  for (name in names(data_lists)) {
    cat("Data for name:", name, "\n")
    cat("Spectral distances calculated:",data_lists[[name]],"\n\n")
    indices <- list(...)
    if (length(indices) == 0) {
      cat("No values provided in the ellipsis.\n")
      return(NULL)
    }
    sum_of_indices <- sum(unlist(indices))
    final_distances <- c()
    for(distance in data_lists[[name]]){
      new_distance <- (distance)/(sum_of_indices)
      final_distances <- c(final_distances,new_distance)
    }
    final_distances_list[[name]] <- final_distances
    cat("Spectral distances normalized:",final_distances,"\n\n")
    
  }
  plot_box_diagrams_normalized("Indices-Criteria",final_distances_list,plot_filename)
}

plot_box_diagrams <- function(plot_filename="",data_lists,network)
{
  spectral_data <- data.frame(
    Model = rep(names(data_lists), times = sapply(data_lists, length)),
    SpectralDistance = unlist(data_lists)
  )
  if(nchar(plot_filename) == 0)
  {
    print(
      ggplot(spectral_data, aes(x = Model, y = SpectralDistance, fill = Model)) +
        geom_boxplot() +
        labs(x = "Model Types", y = "Spectral Distance", title = paste0("Box Plot of Spectral Distances ", network)) +
        scale_fill_discrete(name = "Model Types") +
        theme_minimal())
  }
  else
  {
    plot_filename <- file.path(plot_filename,"not_normalized.png")
    gg<-ggplot(spectral_data, aes(x = Model, y = SpectralDistance, fill = Model)) +
      geom_boxplot() +
      labs(x = "Model Types", y = "Spectral Distance", title = paste0("Box Plot of Spectral Distances ", network)) +
      scale_fill_manual(values = color_palette) +
      theme_minimal()+theme(
        panel.background = element_blank(),   # Remove panel background color
        plot.background = element_rect(fill = "white"),  # Set plot background color
        legend.background = element_rect(fill = "white"),  # Set legend background color
        legend.text = element_text(color = "black"),   # Set legend text color
        legend.title = element_text(color = "black"),  # Set legend title color
        axis.title = element_text(color = "black"),    # Set axis title color
        axis.text = element_text(color = "black"),     # Set axis text color
        strip.text = element_text(color = "black")     # Set strip text color (facet labels)
      )
    ggsave(plot_filename,gg)
  }
  
  
}
plot_box_diagrams_normalized <- function(criteria_name,data_spr_normal,plot_filename="")
{
  spectral_data <- data.frame(
    Model = rep(names(data_spr_normal), times = sapply(data_spr_normal, length)),
    SpectralDistance = unlist(data_spr_normal)
  )
  if(nchar(plot_filename) == 0)
  {
    print(
      ggplot(spectral_data, aes(x = Model, y = SpectralDistance, fill = Model)) +
        geom_boxplot() +
        labs(x = "Model Types", y = "Spectral Distance", title = paste("Box Plot of Spectral Distances - ",criteria_name)) +
        scale_fill_discrete(name = "Model Types") +
        theme_minimal())
  }
  else
  {
      plot_filename <- file.path(plot_filename, paste0(criteria_name,".png"))
      gg <- ggplot(spectral_data, aes(x = Model, y = SpectralDistance, fill = Model)) +
        geom_boxplot() +
        labs(x = "Model Types", y = "Spectral Distance", title = paste("Box Plot of Spectral Distances - ",criteria_name)) +
        scale_fill_manual(values = color_palette) +
        theme_minimal()+theme(
          panel.background = element_blank(),   # Remove panel background color
          plot.background = element_rect(fill = "white"),  # Set plot background color
          legend.background = element_rect(fill = "white"),  # Set legend background color
          legend.text = element_text(color = "black"),   # Set legend text color
          legend.title = element_text(color = "black"),  # Set legend title color
          axis.title = element_text(color = "black"),    # Set axis title color
          axis.text = element_text(color = "black"),     # Set axis text color
          strip.text = element_text(color = "black"),     # Set strip text color (facet labels)
          plot.title = element_text(size=12)
        )
      ggsave(plot_filename,gg)
  }
  
  
}
get_minmax_normalized_values_models <- function(data_lists_all_graphs_x,data_lists_all_graphs_y)
{
  min_value <- Inf
  min_value_laplacian <- Inf
  max_value <- 0
  max_value_laplacian <- 0
  final_distances_list <- list()
  print("######## MIN-MAX NORMALIZATION DATA #############")
  for (network in names(data_lists_all_graphs_x)) {
    min_value <- Inf
    min_value_laplacian <- Inf
    max_value <- 0
    max_value_laplacian <- 0
    for(nullmodel in names(data_lists_all_graphs_x[[network]]))
    {
      vector_min_adjacency <- min(data_lists_all_graphs_x[[network]][[nullmodel]])# Find the minimum value in the current vector
      vector_max_adjacency <- max(data_lists_all_graphs_x[[network]][[nullmodel]])
      vector_min_laplacian <- min(data_lists_all_graphs_y[[network]][[nullmodel]])
      vector_max_laplacian <- max(data_lists_all_graphs_y[[network]][[nullmodel]])
      if (vector_min_adjacency < min_value) {
        min_value <- vector_min_adjacency  # Update the minimum value if a smaller value is found
      }
      if(vector_max_adjacency > max_value)
      {
        max_value <- vector_max_adjacency
      }
      if(vector_min_laplacian < min_value_laplacian)
      {
        min_value_laplacian <- vector_min_laplacian
      }
      if(vector_max_laplacian > max_value_laplacian)
      {
        max_value_laplacian <- vector_max_laplacian
      }
    }
    for(nullmodel in names(data_lists_all_graphs_x[[network]]))
    {
      data_lists_all_graphs_x[[network]][[nullmodel]] <- (data_lists_all_graphs_x[[network]][[nullmodel]] - min_value) / (max_value - min_value)
      data_lists_all_graphs_y[[network]][[nullmodel]] <- (data_lists_all_graphs_y[[network]][[nullmodel]] - min_value_laplacian) / (max_value_laplacian - min_value_laplacian)
      
    }
    
  }
  return(list(x_values=data_lists_all_graphs_x,y_values=data_lists_all_graphs_y))
}
get_covariance_normalized_values <- function(data_lists_all_graphs_x,data_lists_all_graphs_y)
{
  for(network in names(data_lists_all_graphs_x))
  {
    max_length <- max(sapply(data_lists_all_graphs_x[[network]], length))
    max_length_laplacian <- max(sapply(data_lists_all_graphs_y[[network]], length))
    # Fill in missing values with NA to make all vectors of equal length
    filled_variables <- lapply(data_lists_all_graphs_x[[network]], function(x) {
      if (length(x) < max_length) {
        c(x, rep(NA, max_length - length(x)))
      } else {
        x
      }
    })
    filled_variables_laplacian <- lapply(data_lists_all_graphs_y[[network]], function(x) {
      if (length(x) < max_length_laplacian) {
        c(x, rep(NA, max_length_laplacian - length(x)))
      } else {
        x
      }
    })
    
    # Create a data frame from the list
    df <- data.frame(filled_variables)
    df_laplacian <- data.frame(filled_variables_laplacian)
    cov_matrix <- cov(df,use = "pairwise.complete.obs")
    cov_matrix_laplacian <- cov(df_laplacian,use = "pairwise.complete.obs")
    eigenvalues <- eigen(cov_matrix)$values
    eigenvalues_laplacian <- eigen(cov_matrix_laplacian)$values
    sum_eigenvalues <- sum(eigenvalues)
    sum_eigenvalues_laplacian <- sum(eigenvalues_laplacian)
    for (model in names(data_lists_all_graphs_x[[network]])) {
      data_lists_all_graphs_x[[network]][[model]] <- 1 - exp(-data_lists_all_graphs_x[[network]][[model]] / sum_eigenvalues)
      data_lists_all_graphs_y[[network]][[model]] <- 1 - exp(-data_lists_all_graphs_y[[network]][[model]] / sum_eigenvalues_laplacian)
      
    }
  }
  return(list(x_values=data_lists_all_graphs_x,y_values=data_lists_all_graphs_y))
}
plot_clustering_spd_means <- function(data_lists_all_graphs_x,data_lists_all_graphs_y,option="Min_Max")
{
  if(option=="Min_Max")
  {
    result <- get_minmax_normalized_values_models(data_lists_all_graphs_x,data_lists_all_graphs_y)
    string_needer <- "MinMax_Allnetworks"
  }
  else if(option=="Covariance")
  {
    result <- get_covariance_normalized_values(data_lists_all_graphs_x,data_lists_all_graphs_y)
    string_needer <- "Covariance_Allnetworks"
  }
  x_values_normalized <- result$x_values
  y_values_normalized <- result$y_values
  df <- data.frame(
    "Network" = character(0),
    "Model" = character(0),
    "Adjacency_value" = numeric(0),  
    "Laplacian_value" = numeric(0)
  )
  df_not_weighted <- data.frame(
    "Network" = character(0),
    "Model" = character(0),
    "Adjacency_value" = numeric(0),  
    "Laplacian_value" = numeric(0)
  )
  for(network in names(x_values_normalized))
  {
    result_analysis <- analyze_network(paste(network,".csv",sep=""), directory = datadir, guild_a = "Plant", guild_b = "Pollinator", only_NODF = TRUE)
    weighted_network <- sum(result_analysis$matrix > 1)>0
    if(weighted_network)
    {
      for(nullmodel in names(x_values_normalized[[network]]))
      {
        df <- rbind(
          df,
          data.frame(
            "Network" = network,
            "Model" = nullmodel,
            "Adjacency_value" = mean(x_values_normalized[[network]][[nullmodel]]),
            "Laplacian_value" = mean(y_values_normalized[[network]][[nullmodel]])
          )
        )
      }
    }
    else
    {
      for(nullmodel in names(x_values_normalized[[network]]))
      {
        df_not_weighted <- rbind(
          df_not_weighted,
          data.frame(
            "Network" = network,
            "Model" = nullmodel,
            "Adjacency_value" = mean(x_values_normalized[[network]][[nullmodel]]),
            "Laplacian_value" = mean(y_values_normalized[[network]][[nullmodel]])
          )
        )
      }
    }
  }
  point_size <- 1.5
  text_size <- 1.5
    if(option=="Covariance")
      {
        gg <- scatterplot <- ggplot(df, aes(x = `Laplacian_value`, y = `Adjacency_value`, color = `Model`)) +
          geom_point(size = point_size) +
          scale_color_manual(values = color_palette) +
          geom_text_repel(aes(label = `Network`),max.overlaps = Inf,size=text_size) +  # Add network names as text annotations
          labs(
            x = "Laplacian value",
            y = "Adjacency value",
            title = paste0(paste("Scatterplot of means - Weighted Network - ",option)," Criteria")
          )+
          theme(legend.position = "bottom", plot.margin = margin(15, 15, 15, 15))
        folder_cluster_path <- file.path(cluster_path,paste0(string_needer,"weighted.png"))
        ggsave(folder_cluster_path,gg)
        
        gg <- scatterplot <- ggplot(df_not_weighted, aes(x = `Laplacian_value`, y = `Adjacency_value`, color = `Model`)) +
          geom_point(size = point_size) +
          scale_color_manual(values = color_palette) +
          geom_text_repel(aes(label = `Network`),max.overlaps = Inf,size=text_size) +  # Add network names as text annotations
          labs(
            x = "Laplacian value",
            y = "Adjacency value",
            title = paste0(paste("Scatterplot of means - Not-Weighted Network - ",option)," Criteria")
          )+
          theme(legend.position = "bottom", plot.margin = margin(15, 15, 15, 15))
        folder_cluster_path <- file.path(cluster_path,paste0(string_needer,"not-weighted.png"))
        ggsave(folder_cluster_path,gg)
      }
      else if(option=="Min_Max")
      {
        gg <- scatterplot <- ggplot(df, aes(x = `Laplacian_value`, y = `Adjacency_value`, color = `Model`)) +
          geom_point(size = point_size) +
          scale_color_manual(values = color_palette) +
          geom_text_repel(aes(label = `Network`),max.overlaps = Inf,size=text_size) +  # Add network names as text annotations
          labs(
            x = "Laplacian value",
            y = "Adjacency value",
            title = paste0(paste("Scatterplot of means - Weighted Network - ",option)," Criteria")
          )+
          theme(legend.position = "bottom", plot.margin = margin(15, 15, 15, 15))
        folder_cluster_path <- file.path(cluster_path,paste0(string_needer,"weighted.png"))
        ggsave(folder_cluster_path,gg)
        
        gg <- scatterplot <- ggplot(df_not_weighted, aes(x = `Laplacian_value`, y = `Adjacency_value`, color = `Model`)) +
          geom_point(size = point_size) +
          scale_color_manual(values = color_palette) +
          geom_text_repel(aes(label = `Network`),max.overlaps = Inf,size=text_size) +  # Add network names as text annotations
          labs(
            x = "Laplacian value",
            y = "Adjacency value",
            title = paste0(paste("Scatterplot of means - Not-Weighted Network - ",option)," Criteria")
          )+
          theme(legend.position = "bottom", plot.margin = margin(15, 15, 15, 15))
        folder_cluster_path <- file.path(cluster_path,paste0(string_needer,"not-weighted.png"))
        ggsave(folder_cluster_path,gg)
      }
      
   
  
  
}
check_bipartite_connectedness <- function(adjacency_matrix)
{
  graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
  
 
  return(is_bipartite(graph))
}
change_maxmin_values <- function(ad_spd,lapl_spd)#,lapl_nm_spd)
{
  if(ad_spd > maximum_spectral_distance_found)
  {
    maximum_spectral_distance_found <<- ad_spd
  }
  
  if(lapl_spd > maximum_spectral_distance_found_laplacian)
  {
    maximum_spectral_distance_found_laplacian <<- lapl_spd
  }
  
  # if(lapl_nm_spd > maximum_spectral_distance_found_laplacian_nm)
  # {
  #   maximum_spectral_distance_found_laplacian_nm <<- lapl_nm_spd
  # }
  
  if(ad_spd < minimum_spectral_distance_found)
  {
    minimum_spectral_distance_found <<- ad_spd
  }
  if(lapl_spd < minimum_spectral_distance_found_laplacian)
  {
    minimum_spectral_distance_found_laplacian <<- lapl_spd
  }
  # if(lapl_nm_spd < minimum_spectral_distance_found_laplacian_nm)
  #   minimum_spectral_distance_found_laplacian_nm <<- lapl_nm_spd
  
}

covariance_eigenvalues_normalization <- function(plot_filename="", data_lists_x,network)
{
  max_length <- max(sapply(data_lists_x, length))
  
  # Fill in missing values with NA to make all vectors of equal length
  filled_variables <- lapply(data_lists_x, function(x) {
    if (length(x) < max_length) {
      c(x, rep(NA, max_length - length(x)))
    } else {
      x
    }
  })
  
  # Create a data frame from the list
  df <- data.frame(filled_variables)
  cov_matrix <- cov(df,use = "pairwise.complete.obs")
  eigenvalues <- eigen(cov_matrix)$values
  sum_eigenvalues <- sum(eigenvalues)
  final_distances_list <- list()
  print("######## COVARIANCE NORMALIZATION DATA #############")
  for (name in names(data_lists)) {
    cat("Data for name:", name, "\n")
    cat("Spectral distances calculated:",data_lists[[name]],"\n\n")
    final_distances <- c()
    for(distance in data_lists[[name]]){
      new_distance <- 1 - exp(-distance / sum_eigenvalues)
      final_distances <- c(final_distances,new_distance)
      
    }
    final_distances_list[[name]] <- final_distances
    cat("Spectral distances normalized:",final_distances,"\n\n")
  }
  plot_box_diagrams_normalized(paste0(network,"-Covariance-Eigenvalues-Criteria"),final_distances_list,plot_filename)
  
}

plot_means_histogram <- function(data_lists_all_graphs_x,data_lists_all_graphs_y,option="Min_Max")
{
  if(option=="Min_Max")
  {
    result <- get_minmax_normalized_values_models(data_lists_all_graphs_x,data_lists_all_graphs_y)
    string_needer <- "MinMax_Allnetworks"
  }
  else if(option=="Covariance")
  {
    result <- get_covariance_normalized_values(data_lists_all_graphs_x,data_lists_all_graphs_y)
    string_needer <- "Covariance_Allnetworks"
  }
  x_values_normalized <- result$x_values
  y_values_normalized <- result$y_values
  df <- data.frame(
    "Network" = character(0),
    "Model" = character(0),
    "Adjacency_value" = numeric(0),  
    "Laplacian_value" = numeric(0)
  )
  df_not_weighted <- data.frame(
    "Network" = character(0),
    "Model" = character(0),
    "Adjacency_value" = numeric(0),  
    "Laplacian_value" = numeric(0)
  )
  for(network in names(x_values_normalized))
  {
    result_analysis <- analyze_network(paste(network,".csv",sep=""), directory = datadir, guild_a = "Plant", guild_b = "Pollinator", only_NODF = TRUE)
    weighted_network <- sum(result_analysis$matrix > 1)>0
    if(weighted_network)
    {
      for(nullmodel in names(x_values_normalized[[network]]))
      {
        df <- rbind(
          df,
          data.frame(
            "Network" = network,
            "Model" = nullmodel,
            "Adjacency_value" = mean(x_values_normalized[[network]][[nullmodel]]),
            "Laplacian_value" = mean(y_values_normalized[[network]][[nullmodel]])
          )
        )
      }
    }
    else
    {
      for(nullmodel in names(x_values_normalized[[network]]))
      {
        df_not_weighted <- rbind(
          df_not_weighted,
          data.frame(
            "Network" = network,
            "Model" = nullmodel,
            "Adjacency_value" = mean(x_values_normalized[[network]][[nullmodel]]),
            "Laplacian_value" = mean(y_values_normalized[[network]][[nullmodel]])
          )
        )
      }
    }
  }
    lines_divider_size <- 0.1
    sequence_divider_size <- 0.1
    binwidth <- 0.1
    gg <- ggplot(data = df, aes(x = `Adjacency_value`, fill = `Model`)) +
      geom_histogram(binwidth = binwidth, position = "identity",boundary = 0) +  # Adjust the binwidth as needed
      labs(
        title = paste0(paste0("Histogram of Values by Model - Weighted - Adjacency - ",option), " - Criteria"),
        x = "Value",
        y = "Frequency"
      ) +
      scale_fill_manual(values = color_palette)+
      facet_wrap(~Model, ncol=max(1,ceiling(sqrt(length(unique(df_not_weighted$Model))))),scales="free")+
      theme(strip.text = element_text(size = 12)) +  # Adjust strip text size
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_x_continuous(
        breaks = seq(0, 1, by = sequence_divider_size),  # Set custom breaks at 0.05 intervals
        limits = c(0, 1)  # Set the x-axis limits from 0 to 1
      )+
      geom_vline(xintercept = seq(0.05, 0.95, by = lines_divider_size), color = "lightgray", linetype = "dashed")
    
    folder_histogram_path <- file.path(histogram_path,paste0(string_needer,"Weighted.png"))
    ggsave(folder_histogram_path,gg)
    
    gg <- ggplot(data = df, aes(x = `Laplacian_value`, fill = `Model`)) +
      geom_histogram(binwidth = binwidth, position = "identity",boundary = 0) +  # Adjust the binwidth as needed
      labs(
        title = paste0(paste0("Histogram of Values by Model - Weighted - Laplacian - ",option), " - Criteria"),
        x = "Value",
        y = "Frequency"
      ) +
      scale_fill_manual(values = color_palette)+
      facet_wrap(~Model, ncol=max(1,ceiling(sqrt(length(unique(df_not_weighted$Model))))),scales="free")+
      theme(strip.text = element_text(size = 12)) +  # Adjust strip text size
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_x_continuous(
        breaks = seq(0, 1, by = sequence_divider_size),  # Set custom breaks at 0.05 intervals
        limits = c(0, 1)  # Set the x-axis limits from 0 to 1
      )+
      geom_vline(xintercept = seq(0.05, 0.95, by = lines_divider_size), color = "lightgray", linetype = "dashed")
    
    folder_histogram_path <- file.path(histogram_path,paste0(string_needer,"Laplacian_Weighted.png"))
    ggsave(folder_histogram_path,gg)
    
    gg <- ggplot(data = df_not_weighted, aes(x = `Adjacency_value`, fill = `Model`)) +
      geom_histogram(binwidth = binwidth, position = "identity",boundary = 0) +  # Adjust the binwidth as needed
      labs(
        title =  paste0(paste0("Histogram of Values by Model - Not-Weighted - Adjacency - ",option), " - Criteria"),
        x = "Value",
        y = "Frequency"
      ) +
      scale_fill_manual(values = color_palette)+
      facet_wrap(~Model, ncol=max(1,ceiling(sqrt(length(unique(df_not_weighted$Model))))),scales="free_x")+
      theme(strip.text = element_text(size = 12)) +  # Adjust strip text size
      theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
      scale_x_continuous(
        breaks = seq(0, 1, by = sequence_divider_size),  # Set custom breaks at 0.05 intervals
        limits = c(0, 1)  # Set the x-axis limits from 0 to 1
      )+
      geom_vline(xintercept = seq(0.05, 0.95, by =lines_divider_size), color = "lightgray", linetype = "dashed")
    
    folder_histogram_path <- file.path(histogram_path,paste0(string_needer,"Adjacency_Notweighted.png"))
    ggsave(folder_histogram_path,gg)
    
    gg <- ggplot(data = df_not_weighted, aes(x = `Laplacian_value`, fill = `Model`)) +
      geom_histogram(binwidth = binwidth, position = "identity",boundary = 0) +  # Adjust the binwidth as needed
      labs(
        title =  paste0(paste0("Histogram of Values by Model - Not-Weighted - Laplacian - ",option), " - Criteria"),
        x = "Value",
        y = "Frequency"
      ) +
      scale_fill_manual(values = color_palette)+
      facet_wrap(~Model, ncol=max(1,ceiling(sqrt(length(unique(df_not_weighted$Model))))),scales="free")+
      theme(strip.text = element_text(size = 12)) +  # Adjust strip text size
      theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
      scale_x_continuous(
        breaks = seq(0, 1, by = sequence_divider_size),  # Set custom breaks at 0.05 intervals
        limits = c(0, 1)  # Set the x-axis limits from 0 to 1
      )+
      geom_vline(xintercept = seq(0.05, 0.95, by =lines_divider_size), color = "lightgray", linetype = "dashed")
    
    folder_histogram_path <- file.path(histogram_path,paste0(string_needer,"Laplacian_NotWeighted.png"))
    ggsave(folder_histogram_path,gg)

}
################ Main code ##################
color_palette <- c("RND" = "#1f77b4", "SYTR" = "#ff7f0e", "SHUFFLE" = "#2ca02c", 
                   "VAZ" = "#d62728", "SWAP" = "#9467bd", "R2D" = "#8c564b", 
                   "MGEN" = "#e377c2", "DU_05" = "#7f7f7f")
datadir <- "data/"
csv_directory <- "./data"

# Get a list of all .csv files in the directory
csv_files <- list.files(csv_directory, pattern = "\\.csv$", full.names = FALSE)

# Remove the .csv extension and store the file names in a vector

file_names <- sub("\\.csv$", "", csv_files) 

# Print the vector of file names
print(file_names)
## If you want to specify a model, substitute this following line 
## by a vector of models locations in file system as shown below
lnetw <- file_names
  #c("M_PL_001","M_PL_002","M_PL_003","M_PL_004","M_PL_005","M_PL_006","M_PL_007","M_PL_008","M_PL_009","M_PL_010","M_PL_011","M_PL_012","M_PL_013","M_PL_014","M_PL_016")
 
  #c("M_PL_007","M_PL_008","M_PL_009","M_PL_010","M_PL_011","M_PL_012","M_PL_013","M_PL_014","M_PL_016")
lnetw
num_experiments <- 1000
ppi <- 300
INCOMPLETE_MODEL_VALUE <- -1
plottofile <- TRUE
nmagnitudes <- list("spect_rad","adj_energy","lpl_energy")#,"algebraic_connectivity")

maximum_spectral_distance_found <- 0
minimum_spectral_distance_found <- Inf 
maximum_spectral_distance_found_laplacian <- 0
minimum_spectral_distance_found_laplacian <- Inf 
maximum_spectral_distance_found_laplacian_nm <- 0
minimum_spectral_distance_found_laplacian_nm <- Inf 
### Folders to where the plots will be saved
base_path <- "./box_plots_data_per_network"
cluster_path <- "./clustering_data_graph"
histogram_path <- "./histogram_data_graph"


fract_dummy <- c(5)
folder_path_list <- list()

for (netw in (lnetw)){
  
  
  data_lists <<- list()
  data_lists_laplacian <<- list()
  data_lists_laplacian_normalized <<- list()
  spectral_distances_bounds <<- list()
  spectral_distances_bounds_laplacian <<- list()
  spectral_distances_bounds_laplacian_normalized <<- list()
  result_analysis <- analyze_network(paste(netw,".csv",sep=""), directory = datadir, guild_a = "Plant", guild_b = "Pollinator", only_NODF = TRUE)
  nodes_a <- result_analysis$num_guild_a
  nodes_b <- result_analysis$num_guild_b
  num_links <- result_analysis$links
  num_nodes <- nodes_a+nodes_b
  adj_sq_matrix <- sq_adjacency(result_analysis$matrix, nodes_a, nodes_b)
  requirements <- check_bipartite_connectedness(adj_sq_matrix)
  lapl_matrix <- 0-adj_sq_matrix
  degrees <- rowSums(adj_sq_matrix)
  for (i in 1:nrow(lapl_matrix))
    lapl_matrix[i,i] <- degrees[i]
  
  matrix_lapl_norm <- diag(1 / sqrt(degrees)) %*% lapl_matrix %*% diag(1 / sqrt(degrees))
  network_spectrum_laplacian_norm <- eigen(matrix_lapl_norm)
  network_spectrum_laplacian <- eigen(lapl_matrix)
  network_spectrum = eigen(adj_sq_matrix)
  save_parameters_laplacian(length(network_spectrum_laplacian$values),max(network_spectrum_laplacian$values))
  save_parameters_laplacian_normalized(length(network_spectrum_laplacian_norm$values),max(network_spectrum_laplacian_norm$values))
  
  ## Save for bounds normalization criteria Laplacian
  # save_parameters_laplacian(length(network_spectrum$values),max(network_spectrum$values))
  
  weighted_network <- sum(result_analysis$matrix > 1)>0
  if (weighted_network){
    mnames <- c("RND","SYTR","SHUFFLE","VAZ","SWAP")#,"R2D")
    folder_path <- file.path(base_path, "Weighted Networks")
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    folder_path_list[[netw]] <- file.path(folder_path,netw)
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
  }
  else
  {
    mnames <- c("RND","MGEN","SHUFFLE","VAZ",paste0("DU_",sprintf("%02d",fract_dummy)))
    folder_path <- file.path(base_path, "Not-Weighted Networks")
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    folder_path_list[[netw]] <- file.path(folder_path,netw)
    if (!dir.exists(folder_path_list[[netw]])) {
      dir.create(folder_path_list[[netw]])
    } 
  }
    
  model_full <- create_models_list(mnames,NULL)
  found_model_full <- create_models_list(mnames,FALSE)
  #spectral_distance_list <- c()
  nodes_guild_a <- seq(1,result_analysis$num_guild_a)  # guild a top rows/left columns
  nodes_guild_b <- seq(1,result_analysis$num_guild_b)  # guild b bottom rows/right columns
  for (k in 1:num_experiments) {
    for (tmodel in mnames)
    {
      print(tmodel)
      p <- gen_null_model(tmodel,result_analysis,nodes_guild_a,nodes_guild_b,
                          num_links,num_nodes,nmagnitudes)
      
      ## Comprobamos que en la matriz de adyacencia al menos haya una interacción en cada fila y en cada columna
      # - Si existe al menos una columna o fila cuya suma sea 0, no se puede calcular el espectro de la matriz laplaciana normalizada
      # - En cualquier otro caso, con que al menos una de las columnas o filas sea distinta de cero,
      #   se pueden calcular la de adyacencia y la laplaciana
        fmgr <- find_model_good(p$incidmatrix) ##
        if (fmgr$found_model_full){
          null_adjacency <- sq_adjacency(p$incidmatrix,nodes_a,nodes_b) ## ASK if nodes_a, nodes_b are correct
          null_spectrum <- eigen(null_adjacency)
          
          null_matrix_lapl <- 0-null_adjacency
          null_degrees <- rowSums(null_adjacency)
          for (i in 1:nrow(null_matrix_lapl)){
            null_matrix_lapl[i,i] <- null_degrees[i]
          }
          null_spectrum_laplacian <- eigen(null_matrix_lapl)
          
         
          
          spectral_distance <- sqrt(sum((network_spectrum$values - null_spectrum$values)^2))
          spectral_distance_laplacian <- sqrt(sum((network_spectrum_laplacian$values - null_spectrum_laplacian$values)^2))
          
          change_maxmin_values(spectral_distance,spectral_distance_laplacian)
          ## Save for bounds normalization criteria Adjacency
          bounds<-calculate_bounds_adjacency(adj_sq_matrix,null_adjacency)
          #bounds_laplacian <- calculate_bounds_adjacency(lapl_matrix,null_matrix_lapl)
          
          
          calculate_and_save_spdbound(tmodel,bounds,spectral_distance,length(network_spectrum$values))
          #calculate_and_save_spdbound_laplacian(tmodel,bounds_laplacian,spectral_distance_laplacian,length(network_spectrum_laplacian$values))
         
          ## Save for bounds normalization criteria Laplacian Normalized
          
          #null_matrix_lapl_norm <- diag(1 / sqrt(null_degrees)) %*% null_matrix_lapl %*% diag(1 / sqrt(null_degrees))
          #null_matrix_lapl_norm_spectrum <- eigen(null_matrix_lapl_norm)
          #spectral_distance_laplacian_norm <- sqrt(sum((network_spectrum_laplacian_norm$values-null_matrix_lapl_norm_spectrum$values )^2))
          # change_maxmin_values(spectral_distance,spectral_distance_laplacian,spectral_distance_laplacian_norm)
          # bounds_laplacian_normalized <- calculate_bounds_adjacency(matrix_lapl_norm, null_matrix_lapl_norm)
          # calculate_and_save_spdbound_laplacian_normalized(tmodel,bounds_laplacian_normalized,spectral_distance_laplacian_norm,length(network_spectrum_laplacian_norm$values))
          # add_data_laplacian_normalized(tmodel,spectral_distance_laplacian_norm)
          
          
          ## Save spectral distance according to the null model type
          add_data(tmodel,spectral_distance)
          add_data_laplacian(tmodel,spectral_distance_laplacian)
          
          
          
          
        }
      
      
    }
  }
  data_lists_all_graphs <- add_data_all_graphs(netw,data_lists_all_graphs,data_lists)
  data_lists_all_graphs_laplacian <- add_data_all_graphs(netw,data_lists_all_graphs_laplacian,data_lists_laplacian)
  # data_lists_all_graphs_laplacian_normalized <- add_data_all_graphs(netw,data_lists_all_graphs_laplacian_normalized,data_lists_laplacian_normalized)
  laplacian_folder_path <- file.path(folder_path_list[[netw]], "Laplacian")
  adjacency_folder_path <- file.path(folder_path_list[[netw]], "Adjacency")
  # laplacian_normalized_folder_path <- file.path(folder_path, "NM_Laplacian")
  
  
  plot_box_diagrams(adjacency_folder_path,data_lists,netw)
  plot_box_diagrams(laplacian_folder_path,data_lists_laplacian,netw)
  # plot_box_diagrams(laplacian_normalized_folder_path,data_lists_laplacian_normalized)
  
  min_max_normalization(adjacency_folder_path,data_lists,netw)
  min_max_normalization(laplacian_folder_path,data_lists_laplacian,netw)
  # min_max_normalization(laplacian_normalized_folder_path,data_lists_laplacian_normalized)
  
  min_max_normalization_model_independant(adjacency_folder_path,data_lists,netw)
  min_max_normalization_model_independant(laplacian_folder_path,data_lists_laplacian,netw)
  # min_max_normalization_model_independant(laplacian_normalized_folder_path,data_lists_laplacian_normalized)
  
  bounds_normalization(adjacency_folder_path,netw)
  bounds_normalization_laplacian(laplacian_folder_path,network=netw)
  # bounds_normalization_laplacianNormalized(laplacian_normalized_folder_path)
  
  iqr_normalization_all_distances_bynetwork(adjacency_folder_path,data_lists,netw)
  iqr_normalization_all_distances_bynetwork(laplacian_folder_path,data_lists_laplacian,netw)
  # iqr_normalization_all_distances_bynetwork(laplacian_normalized_folder_path,data_lists_laplacian_normalized)
  
  #iqr_normalization(adjacency_folder_path,data_lists)
  #iqr_normalization(laplacian_folder_path,data_lists_laplacian)
  #iqr_normalization(laplacian_normalized_folder_path,data_lists_laplacian_normalized)

  
  covariance_eigenvalues_normalization(adjacency_folder_path,data_lists,netw)
  covariance_eigenvalues_normalization(laplacian_folder_path,data_lists_laplacian,netw)
}

cat("Distancia espectral de adyacencia máxima:" , maximum_spectral_distance_found)
cat("Distancia espectral de adyacencia mínima:" , minimum_spectral_distance_found)
cat("Distancia espectral laplaciana máxima:" , maximum_spectral_distance_found_laplacian)
cat("Distancia espectral laplaciana mínima:" , minimum_spectral_distance_found_laplacian)

plot_clustering_spd_means(data_lists_all_graphs,data_lists_all_graphs_laplacian)
plot_clustering_spd_means(data_lists_all_graphs,data_lists_all_graphs_laplacian,option="Covariance")
plot_means_histogram(data_lists_all_graphs,data_lists_all_graphs_laplacian)
plot_means_histogram(data_lists_all_graphs,data_lists_all_graphs_laplacian,option="Covariance")

#min_max_normalization_all_distances(data_lists_all_graphs,maximum_spectral_distance_found,minimum_spectral_distance_found,"Adjacency")
#min_max_normalization_all_distances(data_lists_all_graphs_laplacian,maximum_spectral_distance_found_laplacian,minimum_spectral_distance_found_laplacian,"Laplacian")

#iqr_normalization_all_distances(data_lists_all_graphs,"Adjacency")
#iqr_normalization_all_distances(data_lists_all_graphs_laplacian,"Laplacian")

#iqr_normalization_all_distances_bymodeltype(data_lists_all_graphs,"Adjacency")
#iqr_normalization_all_distances_bymodeltype(data_lists_all_graphs_laplacian,"Laplacian")



