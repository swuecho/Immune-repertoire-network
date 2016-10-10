#!/usr/bin/ Rscript

library(stringr)
library(reshape)
library(plyr)
library(combinat)

ptm <- proc.time()

source("init_shared.R")



data<-read.csv(all_clusters_csv_file,stringsAsFactors = FALSE)
#remove first column, 
# TODO trt data_without_rowname<-read.csv(all_clusters_csv_file, stringsAsFactors = FALSE, row.names=FALSE)
data<-data[,-1]


subset_list<-unique(data$subset)
####### PART I ########
# get number of clusters by susbset

#TODO: better ways to create a dataframe?
count_cluster_per_subset<-as.data.frame(matrix(nrow=length(subset_list),ncol=2))
colnames(count_cluster_per_subset)<-c('subset','number_of_clusters')

##TODO: use group function to do this --Hao
for (i in 1:length(subset_list)){
  data_sub<-data[data$subset==subset_list[i],]
  count_cluster_per_subset[i,]<-c(subset_list[i],length(unique(data_sub$cluster_membership)))
}

print("number of clusters by subset")
print(count_cluster_per_subset)

count_cluster_per_subset<-count_cluster_per_subset[order(count_cluster_per_subset$subset),]
write.csv(count_cluster_per_subset,file.path(result_dir,paste0(patient, "_count_cluster_per_subset.csv")))



# only keep cluster with node from at least 2 subset
cluster<-unique(data$cluster_membership)
final_clus<-data.frame(subset1=as.character(),subset2=as.character())

for (data_clu in split(data, data$cluster_membership)) {
  cluster_member<-as.character(unique(data_clu$subset))
  if (length(cluster_member)>1){
    sub<-t(data.frame(combn(cluster_member, 2,simplify = T)))
    final_clus<-rbind(final_clus,sub)
  }
  
}


write.table(final_clus,file.path(result_dir,paste0(patient, '_connectivity_rawlist')),sep='\t',row.names=F)

data<-final_clus[order(final_clus$V1,final_clus$V2),]

data2<-data.frame(table(data))
data2<-data2[data2$Freq!=0,]
data2$V1<-as.character(data2$V1)
data2$V2<-as.character(data2$V2)

data2$member<-ifelse(data2$V1<data2$V2,paste(data2$V1,data2$V2,sep=';'),paste(data2$V2,data2$V1,sep=';'))
conn<-ddply(data2[,c(3,4)], .(member), summarize,count = sum(Freq))

write.csv(data2,file.path(result_dir,paste0(patient, '_connectivity_rawlist_summary.csv')),row.names=F)
conn$sourc<-str_replace(str_extract(conn$member,'^.*;'),';','')
conn$target<-str_replace(str_extract(conn$member,';.*$'),';','')


conn<-conn[,-1]
write.table(conn,file.path(result_dir,paste0(patient, '_subset_connectome.txt')),row.names=F,col.names=F,quote =F)

paste(commandArgs(), collapse = " ")
print(proc.time() - ptm)
