#!/usr/bin/ Rscript
library(stringr)
require(plyr)
library(combinat)
library(reshape)

ptm <- proc.time()
source("init_shared.R")

## find the cluster numbers for each subset
print(csf_related_clusters_csv_file)
if (!file.exists(csf_related_clusters_csv_file)) {
  stop("!!!!!!!!!!!!!!!no data for CSF_PB connection files, proablabl there are only CSF data or PB data !!!!!!!!!!!!!!!")
}
## find the cluster numbers for each subset
data<-read.csv(csf_related_clusters_csv_file,stringsAsFactors = FALSE)

if (nrow(data) == 0) {
  stop("!!!!!!!!!!!!!!!no data for CSF!!!!!!!!!!!!!!!")
}

data<-data[,-1]


#CLEANUP this

data$subset<-data$PATIENT
#TODO: remove this? Hao
# build a standard to process sample names, or process sample name in the end
#data$subset<-str_replace(data$subset,paste(patient,'_a_',sep=''),'')
#data$subset<-str_replace(data$subset,'_IGVH','')
#data$subset<-str_replace(data$subset,'_\\d+','')


cluster<-unique(data$cluster_membership)
final_clus<-data.frame(subset1=as.character(),subset2=as.character())
for (data_clu in split(data,data$cluster_membership)) {
  #data_clu<-data[data$cluster_membership==clu,]
  cluster_member<-as.character(unique(data_clu$subset))
  if (length(cluster_member)>1){
    sub<-t(data.frame(combn(cluster_member, 2,simplify = T)))
    final_clus<-rbind(final_clus,sub)
    
  }
  
}




head(final_clus)
data<-final_clus[order(final_clus$V1,final_clus$V2),]
data
str(data)
str(final_clus)
table(final_clus)
data.frame(table(final_clus))

data2<-data.frame(table(data))
data2<-data2[data2$Freq!=0,]
data2
str(data2)
data2$V1<-as.character(data2$V1)
data2$V2<-as.character(data2$V2)

data2$member<-ifelse(data2$V1<data2$V2,paste(data2$V1,data2$V2,sep=';'),paste(data2$V2,data2$V1,sep=';'))
table(data2$member)
head(data2)
conn<-ddply(data2[,c(3,4)], .(member), summarize,count = sum(Freq))

head(conn)
# str(conn)
# data2[data2$member=='CSF-DN.IgG_CSF-SM.IgG',]
# write the connectiveity table to be imported to cytoscape
# write.table(final_clus,paste('Connectivity_LGI1_rawlist',patient, '.txt',sep=''),sep='\t',row.names=F)
# write.table(final_clus,paste('Connectivity_LGI1_rawlist',patient, '.txt',sep=''),sep='\t',row.names=F)
conn$sourc<-str_replace(str_extract(conn$member,'^.*;'),';','')
conn$target<-str_replace(str_extract(conn$member,';.*$'),';','')

head(conn)

conn<-conn[,-1]
write.table(conn,file.path(result_dir,paste0(patient, '_subset_connectome_CSF_related.txt')),row.names=F,col.names=T,quote =F)

paste(commandArgs(), collapse = " ")
print(proc.time() - ptm)
