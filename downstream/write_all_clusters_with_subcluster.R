#!/usr/bin/ Rscript
library(stringr)
library(igraph)
library(rjson)
library(dplyr)

# register multicore for foreach
#library(doMC)
#registerDoMC(cores=10)

ptm <- proc.time()
source("init_shared.R")



el<-read.table(overall_cluster_edge_list_file, sep="\t", header=T,stringsAsFactors = FALSE)
data<-read.table(data_for_clustering, sep="\t", header=T,stringsAsFactors = FALSE)
subsets <- unique(data$subset)
# do quality control here
#rearrange the column, so first column is node_id
vertex<-data[,c("Node_ID","cluster_criteria", "compartment","patient","subset", "vgene","jgene","CDR3","PATIENT","count")]

graph_all<- graph.data.frame(el, directed=F,vertices=vertex)

# Calculate the maximal (weakly or strongly) connected components of a graph
vertex$cluster_membership<-as.factor(clusters(graph_all)$membership)

graph_all<- graph.data.frame(el, directed=F,vertices=vertex)

write.graph(graph_all, all_clusters_gml_file , format = "gml")

write.csv(vertex, all_clusters_csv_file )


# only PB or only CSF(this should be rarely happens)
# no need to get csf related subgraph  
if (length(unique(data$compartment)) < 2 ) {
  stop("!!!!!!!!!!!!!!!no data for CSF in write all cluster.R !!!!!!!!!!!!!!!")
}

# find the clusters with nodes comming from CSF to PB and all CSF nodes
# cluster_list<-as.character(unique(vertex$cluster_membership))
connected_cluster_list<-character()
for (data_cluster_sample in split(vertex,vertex$cluster_membership)){
  # all csf and related node are included. i.e. for CSF, even if there is no connection related to it, 
  # we still included them
  if (length(unique(data_cluster_sample$compartment)) > 1 || unique(data_cluster_sample$compartment)=='CSF') {
    connected_cluster_list<-c(connected_cluster_list,first(data_cluster_sample$cluster_membership))
  }
}

cluster_connecting_CSF_PB<-vertex[vertex$cluster_membership%in%connected_cluster_list,]
#TODO:
# change exported file name
write.csv(cluster_connecting_CSF_PB,csf_related_clusters_csv_file)


subgraph_with_all_egde_list<- graph.data.frame(el, directed=F,vertices=vertex)

##exprot the igraph object into gml file for cytoscape visualization

g_connect<-induced.subgraph(subgraph_with_all_egde_list,which(V(subgraph_with_all_egde_list)$name%in%as.character(cluster_connecting_CSF_PB$Node_ID)))

write.graph(g_connect, csf_related_clusters_gml_file , format = "gml")

#plot the graph
csf_related_plot <- file.path(result_dir, paste0(patient, "_csf_related_graph", ".pdf"))
pdf(csf_related_plot)
print(plot(g_connect))
dev.off()





# cluster which all vertex is from the same subset
# this is too slow
#for (vertex_in_one_cluster in split(vertex,vertex$cluster_membership)){
#   membership <- vertex_in_one_cluster$cluster_membership[1]
#   vertex[vertex$cluster_membership == membership,]$subset_count_in_cluster<- length(unique(vertex_in_one_cluster$subset)) 
#
#}

# membership_with_more_than_one_subset <-c()
# for (vertex_in_one_cluster in split(vertex,vertex$cluster_membership)){
#  
#   membership <- unique(vertex_in_one_cluster$cluster_membership)
#    if(length(unique(vertex_in_one_cluster$subset)) > 1) {
#      append(membership_with_more_than_one_subset, membership)
#    }
# }

# membership_with_more_than_one_subset<-foreach(data_in_subset = split(vertex[1:1000,],vertex$cluster_membership), .combine=c) %dopar% {
#   membership <- unique(vertex_in_one_cluster$cluster_membership)
#   if(length(unique(vertex_in_one_cluster$subset)) > 1) {
#     print(membership)
#     return(membership)
#   }
# }
# find the clusters with nodes comming from CSF to PB and all CSF nodes
# cluster_list<-as.character(unique(vertex$cluster_membership))
membership_with_more_than_one_subset<-character()
for (data_cluster in split(vertex,vertex$cluster_membership)){
  if (length(unique(data_cluster$subset)) > 1) {
    membership_with_more_than_one_subset<-c(membership_with_more_than_one_subset,first(data_cluster$cluster_membership))
  }
}

vertex_remove_innner_subset_cluster<-induced.subgraph(graph_all,
                                                    which(V(graph_all)$cluster_membership %in% membership_with_more_than_one_subset))

cluster_include_at_least_2_subset_gml <- file.path(result_dir, paste0(patient, "_hetero_clusters_", config_suffix, ".gml"))

write.graph(vertex_remove_innner_subset_cluster, cluster_include_at_least_2_subset_gml  , format = "gml")

#### cluster with csf and pb member
membership_with_csf_and_pb <-character()
for (data_cluster in split(vertex,vertex$cluster_membership)){
  if (length(unique(data_cluster$compartment)) > 1) {
    #print(unique(data_cluster$compartment))
    membership_with_csf_and_pb<-c(membership_with_csf_and_pb,first(data_cluster$cluster_membership))
  }
}

vertex_csf_and_pb<-induced.subgraph(graph_all,
                                                      which(V(graph_all)$cluster_membership %in% membership_with_csf_and_pb))

csf_pb_number_of_cluster<-no.clusters(vertex_csf_and_pb)
cluster_include_at_least_csf_and_pb_subset_gml <- file.path(result_dir, paste0(patient, "_hetero_clusters_csf_and_pb_" ,csf_pb_number_of_cluster, '_',config_suffix, ".gml"))

write.graph(vertex_csf_and_pb, cluster_include_at_least_csf_and_pb_subset_gml  , format = "gml")
####

analysis_two_subset <- function (subset_vec,vertex_remove_innner_subset_cluster) {
# find cluster with targeted subsets
graph_with_either_subset<-induced.subgraph(vertex_remove_innner_subset_cluster,
                                                    which(V(vertex_remove_innner_subset_cluster)$subset %in% subset_vec))


# find the clusters with nodes comming from CSF to PB and all CSF nodes
# cluster_list<-as.character(unique(vertex$cluster_membership))

vertex_of_subgraph <- V(graph_with_either_subset)

membership_with_more_than_one_subset<-character()
for (data_cluster in split(vertex_of_subgraph,vertex_of_subgraph$cluster_membership)){
  if (length(unique(data_cluster$subset)) > 1) {
    membership_with_more_than_one_subset<-c(membership_with_more_than_one_subset,first(data_cluster$cluster_membership))
  }
}

graph_between_two_subsets<-induced.subgraph(graph_with_either_subset,
                                                      which(V(graph_with_either_subset)$cluster_membership %in% membership_with_more_than_one_subset))


summary(graph_between_two_subsets)
print(subset_vec)

number_of_clusters <-no.clusters(graph_between_two_subsets)
print(number_of_clusters)

subcluster_dir<- file.path(result_dir, 'subcluster')

if (!file.exists(subcluster_dir )) {
  dir.create(subcluster_dir, recursive=TRUE)
   # keep a copy of config file
}

if(number_of_clusters >0){
  combined_name <- file.path(subcluster_dir, paste(paste(subset_vec[1], subset_vec[2],sep="_x_"), config_suffix,  sep="_"))

    
subset_cluster_gml <- paste0(combined_name, ".gml")
subset_cluster_pdf <- paste0(combined_name, ".pdf")
write.graph(graph_between_two_subsets, subset_cluster_gml , format = "gml")

pdf(subset_cluster_pdf)
print(plot(graph_between_two_subsets))
dev.off()
  
return(number_of_clusters)

}else {
  return(0)
}

}


number_of_subsets <- length(subsets)
mm <- matrix(0,number_of_subsets,number_of_subsets)
rownames(mm)<-subsets
colnames(mm)<-subsets
all_combn<-utils::combn(subsets, m = 2)
for(i in 1:ncol(all_combn)) {
  print(all_combn[,i])
  subset_vec<-all_combn[,i]
  cluster_number <-analysis_two_subset( subset_vec,vertex_remove_innner_subset_cluster)
  print(cluster_number)
  mm[subset_vec[1],subset_vec[2]]<-cluster_number
  mm[subset_vec[2],subset_vec[1]]<-cluster_number 
  print(i)
  print("___-----------___")
}

print(mm)
write.csv(mm,file.path(result_dir,paste0(patient, '_cluster_connection_matrix.csv')),quote =F)

paste(commandArgs(), collapse = " ")
print(proc.time() - ptm)
