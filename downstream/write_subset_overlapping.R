#!/usr/bin/ Rscript

ptm <- proc.time()
source("init_shared.R")

args1 <- file.path(result_dir,paste0(patient,"_count_cluster_per_subset.csv"))
args2 <- file.path(result_dir,paste0(patient, "_connectivity_rawlist_summary.csv"))

if (file.exists(args1) &  file.exists(args2)) {
data1<-read.csv(args1,header=T,stringsAsFactors = FALSE)
data2<-read.csv(args2,header=T,stringsAsFactors = FALSE)

head(data1)
data1<-data1[,-1]

colnames(data2)[1:3]<-c('subset1','subset2','number_overlapping_clusters')
data2<-data2[,-4]
data2$total_number_of_clusters_subset1<-rep(0)
data2$total_number_of_clusters_subset2<-rep(0)
for (i in 1:(dim(data2))[1]){
  data2$total_number_of_clusters_subset1[i]<-data1[data1$subset==data2[i,1],2]
  data2$total_number_of_clusters_subset2[i]<-data1[data1$subset==data2[i,2],2]
}

write.csv(data2,file.path(result_dir,paste0(patient, '_overlapping_count_between_subset.csv')))
}

paste(commandArgs(), collapse = " ")
print(proc.time() - ptm)
print(Sys.time())
