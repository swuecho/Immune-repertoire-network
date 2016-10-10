#!/usr/bin/ Rscript

library(stringdist)
library(stringr)
library(igraph)
# register multicore for foreach
library(doMC)
registerDoMC(cores=12)
ptm <- proc.time()
source("init_shared.R")



data<-read.table(data_for_clustering, sep="\t", header=T,stringsAsFactors=F)

# cluster_criteria is in class charater already
# data$cluster_criteria<-as.character(data$cluster_criteria) # why this?


# find the clusters for each germline usage and build the edge list out of it,
# then combine all the germline together
cluster_by_distance<-function(data, distance_method, distance_cutoff){
  #CDR3<-as.character(data$CDR3)
  CDR3<-data$CDR3
  data_matrix<-stringdistmatrix(CDR3,CDR3,method= distance_method)
  
  # this is not necessary.  --Hao
  # data_matrix<-as.matrix(data_matrix) # when as.matrix is necessay? how to decide if a var is matrix?
  # TODO: distance cutoff, distance method
  
  # two cdr3 is connected if the distance less or equal to 1
  data_matrix<-ifelse(data_matrix<=distance_cutoff,1,0)

  ## remove self connected edge
  for (i in 1:nrow(data_matrix)){
    data_matrix[i,i]<-0
  }
  
  rownames(data_matrix)<-data$Node_ID
  colnames(data_matrix)<-data$Node_ID

  
  ## format the distance matrix into the adjacency matrix 
  data_adj<-graph.adjacency(data_matrix,mode=c("undirected"))
  
  ## get the edge list from the adjacency matrix
  
  data.edgelist<-as.data.frame(get.edgelist(data_adj))
  
  return(data.edgelist)
}

hamming_edge_list_final <-foreach(data_of_cluster= split(data,data$cluster_criteria), .combine='rbind') %dopar% {
  cluster_by_distance(data_of_cluster,config$method, distance_cutoff)
}
# ? sort list?
write.table(hamming_edge_list_final, overall_cluster_edge_list_file ,sep="\t")

# debug information
paste(commandArgs(), collapse = " ")
print(proc.time() - ptm)
