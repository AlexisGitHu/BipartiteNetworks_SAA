library(kcorebip)
library(ggplot2)
library(patchwork)
library(forcats)
source("synthrade_nullmodel.R")

plot_distr_null <- function(df,nvalue,networkname="",title=""){
  plot <- ggplot(data=df)+geom_histogram(aes(x=vals,fill=type),bins=30,position="identity",alpha=0.5)+
    #plot <- ggplot(data=df)+geom_density(aes(x=vals,fill=type),color="transparent",alpha=0.3,bw = "ucv")+
    scale_y_continuous(expand=c(0,0))+xlab("")+theme_bw()+ 
    theme(legend.position = "bottom",legend.text = element_text(size = 8),
          legend.key.size = unit(0.3, 'cm'))+
    guides(fill=guide_legend(nrow=2, byrow=TRUE, title= "Null model"))
  plot <- plot+geom_vline(data = data.frame("val"=nvalue), aes(xintercept = val), 
                          color = "blue", size=0.5,alpha=0.5)+ggtitle(title)
  return(plot)
}

plot_spectr <- function(df,title=""){
  df$type <- fct_relevel(df$type, "RND", "VAZ", "NETWORK")
  if (nrow(df)<20)
    psize <- 2
  else
    psize <- 1
  plot <- ggplot(data=df)+geom_point(aes(x=ind,y=val,color=type),size=psize,
                                     position="identity",alpha=0.5)+
    xlab("")+ylab("Eigen values")+theme_bw()+ggtitle(title)+
    theme_bw()+ theme(legend.position = "bottom",legend.text = element_text(size = 8),
                      legend.key.size = unit(0.3, 'cm'),
                      axis.title.x = element_text(size = 1))+
    guides(color=guide_legend(title= "Null model"))
  return(plot)
}

plot_bipartite <- function(bg, aspect_ratio = 9/35, vframecolor = "grey70", vlabelcex = 4,
                           vsize = 4, vcolor = c("lightblue","pink2"), labelcolor = c("blue","red"),
                           framedisp = FALSE,  color_link = "grey50", vertical = FALSE)
{
  l <- layout.bipartite(bg)
  if (vertical)
    la <- l[, c(2,1)]
  else
    la <- l
  plot.igraph(bg, layout= la,asp=aspect_ratio,vertex.frame.color=vframecolor,
              vertex.label.cex=vlabelcex,vertex.label.color="black",
              vertex.size=vsize, edge.color= color_link, frame = framedisp,
              vertex.color=vcolor[V(bg)$type+1])
}

plot_null_model <- function(modeldata,modeltype,netname,dirn,plotdir,ptfile){
  filenull <- paste0(gsub(".csv","",netname),"_NULL_",modeltype,".csv")
  print(paste("FILENULL",filenull))
  write.csv(modeldata,paste0(dirnulls,filenull))
  result_analysis <- analyze_network(directory = dirn, filenull, guild_a = "Plant", guild_b = "Pollinator", only_NODF = TRUE)
  pgr_null <- ziggurat_graph(dirn,filenull,plotsdir=plotdir,print_to_file = ptfile,show_title = FALSE,weighted_links = "log10")
}

plot_distributions <- function(mresults,magnitude,stitle,vtitle,nname,fdummy){
  mdls <- names(modresults)
  col <- which(names(modresults[[1]])==magnitude)
  df_nulls <- data.frame("vals"=mresults[[1]][,col],"type"=mdls[[1]])
  for (i in 2:length(mdls))
    df_nulls <- rbind(df_nulls,data.frame("vals"= mresults[[mdls[i]]][,col],"type"=mdls[i]))
  pimage <- plot_distr_null(df_nulls,vtitle,networkname=nname,title=sprintf("%s %.2f",stitle,vtitle))
  calc_values <- list("df_nulls"= df_nulls, "pimage" = pimage)
  return(calc_values)
}

LaplEnergy <- function(eigenvalues,num_links,num_nodes){
  return(sum(abs(eigenvalues-2*num_links/num_nodes)))
}

AdjEnergy <- function(eigenvalues){
  return(sum(abs(eigenvalues)))
}

rndnullmodel_build <- function(adjmatrix,nodes_guild_a,nodes_guild_b,num_links,allconnected = FALSE)
{
  mmatrix <- adjmatrix * 0
  # Throw links at random, some nodes may remain disconnected
  if (!allconnected)
    mmatrix[sample(1:length(adjmatrix),num_links,replace = FALSE)] <- 1
  else{
    for (i in 1:nrow(mmatrix))
      mmatrix[i,sample(1:ncol(mmatrix),1)] <- 1   # One link per row
    for (i in 1:ncol(mmatrix))
      mmatrix[sample(1:nrow(mmatrix),1),i] <- 1   # One link per column
    mmatrix[sample(1:length(adjmatrix),num_links-sum(mmatrix),replace = FALSE)] <- 1
  }
  return(mmatrix)
}

# Intercambia un porcentaje de enlaces

dummynullmodel_build <- function(nodes_guild_a,nodes_guild_b,num_links,re,swap_perc)
{
  mmatrix <- re$matrix
  swap_links <- max(1,round(0.01 * swap_perc * num_links))
  ones <- which(result_analysis$matrix > 0)
  zeroes <- which(result_analysis$matrix == 0)
  j = 1
  found <- FALSE
  while ((!found) & (j<10)){
    if (j>1)
      print(paste("dummy",j))
    pmatrix <- mmatrix
    onestozeroes <- sample(ones,swap_links,replace=FALSE)
    zeroestoones <- sample(zeroes,swap_links,replace=FALSE)
    pmatrix[onestozeroes] <- 0
    pmatrix[zeroestoones] <- 1
    j <- j +1
    found <- ((sum(rowSums(pmatrix))>0) & (sum(colSums(pmatrix))>0))
  } 
  return(pmatrix)
}

null_model_process <- function(admatrix,re)
{
  model_matrix <- sq_adjacency(admatrix,re$num_guild_a,re$num_guild_b)
  nullNMP_spect = eigen(model_matrix)
  lapl_nullNMP <- 0-model_matrix
  degreesnullNMP <- rowSums(model_matrix)
  for (i in 1:nrow(lapl_nullNMP))
    lapl_nullNMP[i,i] <- degreesnullNMP[i]
  lpl_spect_nullNMP = eigen(lapl_nullNMP)
  spect_rad_nulls_NMP <- nullNMP_spect$values[1]
  adj_energy_nulls_NMP <- AdjEnergy(nullNMP_spect$values)
  lpl_energy_nulls_NMP <- LaplEnergy(lpl_spect_nullNMP$values,num_links,num_nodes)
  algebraic_connectivity_nulls_NMP <- lpl_spect_nullNMP$values[length(lpl_spect_nullNMP$values)-1]
  
  calc_values <- list("model_matrix" = model_matrix,
                      "null_spect" = nullNMP_spect,
                      "lpl_spect_nulls" = lpl_spect_nullNMP,
                      "spect_rad_nulls" = spect_rad_nulls_NMP,
                      "adj_energy_nulls" = adj_energy_nulls_NMP,
                      "lpl_energy_nulls" = lpl_energy_nulls_NMP,
                      "algebraic_connectivity_nulls" = algebraic_connectivity_nulls_NMP)
  return(calc_values)
}

sq_adjacency <- function(rmatrix, num_guild_a, num_guild_b)
{
  z_matrix <- matrix(0, ncol = num_guild_b, nrow = num_guild_b)
  incidmatrix <- cbind(rmatrix,z_matrix)
  colnames(incidmatrix)<-NULL
  rownames(incidmatrix)<-NULL
  z_matrix <- matrix(0, ncol = num_guild_a, nrow = num_guild_a)
  incidmatrix <- rbind(cbind(z_matrix,t(incidmatrix[,1:num_guild_a])),incidmatrix)
  if (max(incidmatrix)>1){
    # Weighted matrix, binarization
    incidmatrix[incidmatrix>1]=1
  }
  return(incidmatrix)
}

store_model_results <- function(dp,num_links,num_nodes,magnitudes){
  mr <- replicate(length(magnitudes),c())
  names(mr) <- magnitudes
  mr$spect_rad <- dp$null_spect$values[1]
  mr$adj_energy <- AdjEnergy(dp$null_spect$values)
  mr$lpl_energy <- LaplEnergy(dp$lpl_spect_nulls$values,num_links,num_nodes)
  mr$algebraic_connectivity <- dp$lpl_spect_nulls$values[length(dp$lpl_spect_nulls$values)-1]
  return(as.data.frame(mr))
}

gen_null_model <- function(typemodel,resanalysis,nodesga,nodesgb,nlinks,nnodes,magnitudes)
{
  if (typemodel == "RND")
    incidmatrix <- rndnullmodel_build(resanalysis$matrix,nodesga,nodesgb,nlinks)
  if (typemodel == "RNC")
    incidmatrix <- rndnullmodel_build(resanalysis$matrix,nodesga,nodesgb,nlinks,allconnected = TRUE)
  if (typemodel == "VAZ")
    incidmatrix <- nullmodel(resanalysis$matrix,N=1,method="vaz")[[1]]
  if (typemodel == "SWAP")
    incidmatrix <- nullmodel(resanalysis$matrix,N=1,method="swap.web")[[1]]
  if (typemodel == "R2D")
    incidmatrix <- nullmodel(resanalysis$matrix,N=1,method="r2dtable")[[1]]
  if (typemodel == "MGEN")
    incidmatrix <- nullmodel(resanalysis$matrix,N=1,method="mgen")[[1]]
  if (typemodel == "SHUFFLE")
    incidmatrix <- nullmodel(resanalysis$matrix,N=1,method="shuffle.web")[[1]]
  if(typemodel == "SYTR")
    incidmatrix <- SynthTradeNull(resanalysis)$matrix_synth
  
  if (grepl('^DU_', typemodel)){
    fdummy <- as.integer(gsub("DU_","",typemodel))
    incidmatrix <- dummynullmodel_build(nodesga,nodesgb,nlinks,resanalysis,fdummy)
  }
  resp <- null_model_process(incidmatrix,resanalysis)
  return(list("mres"=store_model_results(resp,nlinks,nnodes,magnitudes),"resp"=resp,"incidmatrix"=incidmatrix))
}

find_model_good <- function(matrix){
  # if ((sum(colSums(matrix)==0)==0) &&
  #    (sum(rowSums(matrix)==0)==0))
  if ( (all(sum(colSums(matrix)==0))) &&
       (all(sum(rowSums(matrix)==0))) )
  {
    model_full <- NULL
    found_model_full <- FALSE
  }  else {
    model_full <- matrix
    found_model_full <- TRUE
  } 
  calc_values=list("model_full"=model_full,"found_model_full"=found_model_full)
  return(calc_values)
}
find_model_good2 <- function(matrix){
  row_sums <- rowSums(matrix)
  col_sums <- colSums(matrix)
  # if ((sum(colSums(matrix)==0)==0) &&
  #    (sum(rowSums(matrix)==0)==0))
  if ( (all(row_sums!=0)) &&
       (all(col_sums!=0)) )
  {
    model_full <- matrix
    found_model_full <- TRUE
    
  }  else {
    model_full <- NULL
    found_model_full <- FALSE
  } 
  calc_values=list("model_full"=model_full,"found_model_full"=found_model_full)
  return(calc_values)
}

create_models_list <- function(lnames,init_value)
{
  mylist <- replicate(length(lnames),init_value)
  names(mylist) <- lnames
  return(mylist)
}

store_spect_values <- function(val,typem){
  return(data.frame("val"=val,"ind"=seq(1:length(val)),"type" = typem))
}
