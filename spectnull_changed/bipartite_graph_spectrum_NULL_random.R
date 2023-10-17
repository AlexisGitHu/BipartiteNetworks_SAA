rm(list = ls())
library(kcorebip)
library(ggplot2)
library(patchwork)
library(forcats)
source("synthrade_nullmodel.R")

# In the original files guild a is stored by columns and guild_b by rows

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

create_models_list <- function(lnames,init_value)
{
  mylist <- replicate(length(lnames),init_value)
  names(mylist) <- lnames
  return(mylist)
}

store_spect_values <- function(val,typem){
  return(data.frame("val"=val,"ind"=seq(1:length(val)),"type" = typem))
}

## Escritura de csv


# options( warn = -1 )
datadir <- "data/"

lnetw <- list()
for (k in seq(28,30))
  lnetw <- c(lnetw,sprintf("M_SD_0%02d.csv",k))

filenames <- Sys.glob(paste0(datadir,"*HP*.csv"))
lnetw <- gsub(datadir,"",filenames)

lnetw <- c("RA_HP_006.csv")
lnetw
num_experiments <- 1000
ppi <- 300
INCOMPLETE_MODEL_VALUE <- -1
plottofile <- TRUE
plotzigs <- TRUE
odir <- "plots/"
rdir <- "results/"
dirnulls <- "nullmatrix/"
NetworkMagsFile <- "NetworkMagnitudes.csv"
dir.create(odir, showWarnings = FALSE)
dir.create(rdir, showWarnings = FALSE)
dir.create(dirnulls, showWarnings = FALSE)
fract_dummy <- c(5)
nmagnitudes <- list("spect_rad","adj_energy","lpl_energy")#,"algebraic_connectivity")
for(netw in (lnetw)){
  print("hola")
}
for (netw in (lnetw)){
  # csv_Prueba <-  read.csv('data/RA_HP_047.csv')
  # View(csv_Prueba)
  result_analysis <- analyze_network("M_PL_007.csv", directory = datadir, guild_a = "Plant", guild_b = "Pollinator", only_NODF = TRUE)
  weighted_network <- sum(result_analysis$matrix > 1)>0
  nodes_a <- result_analysis$num_guild_a
  nodes_b <- result_analysis$num_guild_b
  adj_sq_matrix <- sq_adjacency(result_analysis$matrix, nodes_a, nodes_b)
  num_links <- result_analysis$links
  num_nodes <- nodes_a+nodes_b
  adj_spect = eigen(adj_sq_matrix)
  lapl_matrix <- 0-adj_sq_matrix
  adj_energy <- AdjEnergy(adj_spect$values)
  print(sprintf("Sum adjacency spectrum^2 %.2f Energy %.2f",sum((adj_spect$values)^2),adj_energy))
  degrees <- rowSums(adj_sq_matrix)
  for (i in 1:nrow(lapl_matrix))
    lapl_matrix[i,i] <- degrees[i]
  lpl_spect = eigen(lapl_matrix)
  lpl_energy <- LaplEnergy(lpl_spect$values,num_links,num_nodes)
  print(sprintf("Sum Laplacian spectrum %.2f",sum(lpl_spect$values)))
  networkspect <- list("spect_rad"=adj_spect$values[1],"adj_energy"=adj_energy,
                       "lpl_energy"=lpl_energy,
                       "algebraic_connectivity"=lpl_spect$values[length(lpl_spect$values)-1])
  
  if (weighted_network)
    mnames <- c("RND","SYTR","SHUFFLE","VAZ","SWAP")#,"R2D")
  else
    mnames <- c("RND","MGEN","SHUFFLE","VAZ",paste0("DU_",sprintf("%02d",fract_dummy)))
  nmodels <- length(mnames)
  model_full <- create_models_list(mnames,NULL)
  found_model_full <- create_models_list(mnames,FALSE)
  nname <- gsub(".csv","",netw)
  print(nname)
  nullsinfo <- data.frame("spect_rad"=replicate(num_experiments,INCOMPLETE_MODEL_VALUE))
  nullsinfo$adj_energy <- INCOMPLETE_MODEL_VALUE
  nullsinfo$lpl_energy <- INCOMPLETE_MODEL_VALUE
  nullsinfo$algebraic_connectivity <- INCOMPLETE_MODEL_VALUE
  modresults <- lapply(1:nmodels, function(x) nullsinfo)
  names(modresults) <- mnames
  specresults <- lapply(1:nmodels, function(x) c())
  names(specresults) <- mnames
  incidmatrix <- specresults
  datamod <- stack(networkspect)
  datamod$MODEL <- "NETWORK"

  nodes_guild_a <- seq(1,result_analysis$num_guild_a)  # guild a top rows/left columns
  nodes_guild_b <- seq(1,result_analysis$num_guild_b)  # guild b bottom rows/right columns
  for (k in 1:num_experiments) {
    for (tmodel in mnames)
    {
      p <- gen_null_model(tmodel,result_analysis,nodes_guild_a,nodes_guild_B,
                          num_links,num_nodes,nmagnitudes)
      modresults[[tmodel]][k,]<- p$mres
      specresults[[tmodel]] <- p$res
      incidmatrix[[tmodel]] <- p$incidmatrix
      
      if(!found_model_full[[tmodel]]){
        fmgr <- find_model_good(p$incidmatrix)
        if (fmgr$found_model_full){
          model_full[[tmodel]] <- fmgr$model_full
          found_model_full[[tmodel]] <- TRUE
        }
      }
    }
  }
  
  print("BIPARTITE_FOR 1")
  # Save model results
  for (m in names(modresults)){
    sm <- stack(modresults[[m]])
    sm$MODEL <- m
    datamod <- rbind(datamod,sm)
    # if (mark_incomplete_models)
    #   datamod <- datamod[datamod$values != INCOMPLETE_MODEL_VALUE,]
    
  }
  print("BIPARTITE_FOR 2")
  lheader <- length(networkspect)
  numvals <- (nrow(datamod)-lheader)/num_experiments
 
  ntry <- c()
  for (i in 1:lheader)
    ntry <- c(ntry,1:num_experiments)
  dtry <- replicate(lheader,0)
  for (i in 1:length(names(modresults)))
    dtry <- c(dtry,ntry)
  datamod$TRY <- dtry
  write.csv(datamod,paste0(rdir,"MODS_",netw),row.names = FALSE)
  
  testvalues = data.frame("MODEL"=c(),"wstatistic"=c(),"magnitude"=c(),"networkvalue"=c(),
                          "nullavlue"=c(),"pvalue"=c(),"quantile"=c(),"distance"=c(),"reldist"=c())
  for (modeln in mnames)
    for (magnitude in nmagnitudes)
    {
    datos <- datamod[(datamod$MODEL==modeln) & (datamod$ind==magnitude),]$values
    mvalue <- mean(datos)
    datored <- datamod[(datamod$MODEL=="NETWORK") & (datamod$ind==magnitude),]$values
    w <- wilcox.test(x = datos, mu = datored, 
                     alternative = "two.sided") 
    percentile <- ecdf(datos)
    testvalues <- rbind(testvalues,data.frame("MODEL"=modeln,"wstatistic"=w$statistic,
                                             "magnitude"=magnitude, "networkvalue"=datored,
                                             "modelmean"=mvalue, "pvalue"=w$p.value,
                                             "quantile"=100*percentile(datored),
                                             "distance"=datored-mvalue,"reldist"=(datored-mvalue)/datored))
    }
  write.csv(testvalues,paste0(rdir,"TESTVALUES_",netw),row.names = FALSE)
  # Write Network magnitudes to global results file
  network_data <- data.frame("Network"=nname,"NodesA"=nodes_a,"NodesB"=nodes_b,"Links"=num_links,
                 "Weighted"=weighted_network,"Model"="NETWORK")
  network_values <- cbind(network_data,as.data.frame(networkspect))
  for (modelname in mnames){
    network_null <- cbind(network_data,as.data.frame(lapply(modresults[[modelname]],mean)) )
    network_null$Model <- modelname
    network_values <- rbind(network_values,network_null)
  }
  NFile <- paste0(rdir,NetworkMagsFile)
  if (file.exists(NFile)){
      NMags <- read.csv(NFile)
      NMags <- NMags[NMags$Network!=nname,]
      NMags <- rbind(NMags,network_values)
  } 
  else       
      NMags <- network_values
  write.csv(NMags,NFile,row.names = FALSE)
  if (plottofile)
  {
    if (weighted_network)
      mspec <- c("NETWORK","VAZ","SYTR","RND")
    else
      mspec <- c("NETWORK","VAZ","SHUFFLE","RND")
    df_adj_model <- lapply(1:length(mspec), function(x) c())
    names(df_adj_model) <- mspec
    df_lpl_model <- df_adj_model
    
    df_adj_model[["NETWORK"]] <- store_spect_values(adj_spect$values,"NETWORK")
    df_lpl_model[["NETWORK"]] <- store_spect_values(lpl_spect$values,"NETWORK")
    
    for (m in mspec[2:length(mspec)])
    {
      df_adj_model[[m]] <- store_spect_values(specresults[[m]]$null_spect$values,m)
      df_lpl_model[[m]] <- store_spect_values(specresults[[m]]$lpl_spect_null$values,m)
    }
    if (plotzigs)
      ziggurat_graph(datadir,netw,plotsdir=odir,print_to_file = plottofile,show_title = FALSE,weighted_links = "log10")
    pl_adj_spectrum <- plot_spectr(do.call(rbind,df_adj_model),title="Adjacency spectrum")
    pl_lpl_spectrum <- plot_spectr(do.call(rbind,df_lpl_model),title="Laplacian spectrum")
    for (tmodel in mnames)
      if (found_model_full[[tmodel]])
        if (plotzigs)
          plot_null_model(model_full[[tmodel]],tmodel,netw,dirnulls,odir,plottofile)
    print("POST NULL MODEL 1")
    pd_spect_rad <- plot_distributions(modresults,"spect_rad","Spectral radius",adj_spect$values[1],netw,fract_dummy)
    pd_energy <- plot_distributions(modresults,"adj_energy","Energy",adj_energy,netw,fract_dummy)
    pd_lpl_energy <- plot_distributions(modresults,"lpl_energy","Laplacian energy",lpl_energy,netw,fract_dummy)
    wsup <- (pl_adj_spectrum | pl_lpl_spectrum )
    winf <- (pd_spect_rad$pimage | pd_energy$pimage | pd_lpl_energy$pimage)
    wtot <- (wsup / winf) + plot_layout(heights = c(0.5,0.5))
             plot_annotation(title = nname,
                      theme = theme(plot.title = element_text(size = 16,hjust=0.5))) 
    sweight <- "WEIGHTED"
    if (!weighted_network)
      sweight <- "BINARY"
    nfile <- paste0(odir,paste0(nname,"_ALLPLOTS_",sweight))
    png(paste0(nfile,".png"),width=11*ppi,height=7*ppi,res=ppi)
    print(wtot)
    dev.off()
  }
  
}

