#!/usr/bin/ Rscript
### SQL command to collect the patient information
library(rjson)     
library(stringr)
library(dplyr)
library(foreach)

# register multicore for foreach
library(doMC)
registerDoMC(cores=12)
print(Sys.time())
ptm <- proc.time()

############ 
# find_identical_CDR3_across_VJ
# from Shengzhi Wang <wszcas@gmail.com>
# basically what it does it to find out for a particular CDR3 AA sequence with a particular B cell subset, 
# whether there is more than one VJ gene pair assigned to it. 
# Theoritically, it should not because the chance is very low due to the diversity of B cell.  
# If yes, keep the most abundant VJ and remove the rest.  Finally we get a filtered dataset from the original one.
# 
find_identical_CDR3_across_VJ<-function(data_in_subset){
  #cdr3_list<-as.character(unique(data_in_subset$CDR3))
  data_final<-foreach(data_in_subset_with_same_cdr3 = split(data_in_subset,data_in_subset$CDR3), .combine='rbind') %do% {
    #data_in_subset_with_same_cdr3<-data_in_subset[data_in_subset$CDR3==cdr3,]
    n_distinct_VJ<-n_distinct(data_in_subset_with_same_cdr3$cluster_criteria)
    if (n_distinct_VJ > 1) {
      max_read_count<-max(data_in_subset_with_same_cdr3$count)
      #print("multiple VJ combination share same cdr3")
      #print(data_in_subset_with_same_cdr3)
      data_majority<-data_in_subset_with_same_cdr3[data_in_subset_with_same_cdr3$count==max_read_count,]
     
      # the first cluster_criteria has max read_count will be used
      # TODO: the frequency of VJ combination used?
      # 1.
      # if the count of multiple VJ combination are equal to max read cont, 
      # should use the most likely VJ combination
      # 2. 
      # or use the sum of all count?
      if (n_distinct(data_majority$cluster_criteria)>1){
        data_majority<-data_majority[data_majority$cluster_criteria==first(data_majority$cluster_criteria),]
      }
     
    }
    else if (n_distinct_VJ==1) {
      data_majority<-data_in_subset_with_same_cdr3
    }
    
    return (data_majority)
  }
  return(data_final)
}
#################


source("init_shared.R")
# create target dir
if (!file.exists(result_dir)) {
  dir.create(result_dir, recursive=TRUE)
}
data<-read.csv(raw_csv_file, header = TRUE,stringsAsFactors=FALSE)


data <-data[,c("vgene", "jgene", "cdr3", "barcode", "count")]
colnames(data)<-c("vgene","jgene","CDR3","PATIENT","count")  

raw_summary_of_sequence_cout <- data %>% group_by(PATIENT) %>% summarise(count = sum(count)) %>% collect()
raw_summary_of_sequence_cout_path <- file.path(result_dir,paste0(patient, '_sequence_count_summary_raw.txt'))
write.csv(raw_summary_of_sequence_cout, raw_summary_of_sequence_cout_path,row.names = F)

# clean date
# filter
data <- data[data$count>=count_cutoff,]
data <- data[str_count(data$vgene,'IGHV')==1 & str_count(data$jgene,'IGHJ') ==1,]

# merge row with same vgene, jgene, cdr3, and subset
data<- (data %>% group_by(vgene, jgene, CDR3, PATIENT) %>% summarise(count = sum(count)) %>% collect())

data$cluster_criteria<-paste(data$vgene,data$jgene,sep="_")
data$Node_ID<-paste(data$PATIENT,rownames(data),sep="_")

# decide the compartment and subset from barcode

# > unlist(strsplit('31612_b_UNS_1000000_IGG_IGVH','_'))
# [1] "31612"   "b"       "UNS"     "1000000" "IGG"     "IGVH"  
#  57414_a_CSF-DN_15_IGM_IGVH
barcode_splited <- strsplit(data$PATIENT,'_')

# pati
patient_id_vec <- sapply(barcode_splited, '[', 1)
compartment_vec <- sapply(barcode_splited, '[', 3)



# careful, the order of the following three statement matters
# check whole smaple name for CSF 
data$compartment<-ifelse(substr(compartment_vec,1,3) == 'CSF','CSF','PB')
data$patient = patient_id_vec
data$subset<- data$PATIENT


data_clean <-foreach(data_in_subset = split(data,data$subset), .combine='rbind') %dopar% {
  find_identical_CDR3_across_VJ(data_in_subset)
}

clean_summary_of_sequence_cout <- data_clean %>% group_by(PATIENT) %>% summarise(count = sum(count)) %>% collect()
clean_summary_of_sequence_cout_path <- file.path(result_dir,paste0(patient, '_sequence_count_summary_clean.txt'))
write.csv(clean_summary_of_sequence_cout, clean_summary_of_sequence_cout_path,row.names = F)
print(clean_summary_of_sequence_cout)
write.table(data_clean, data_for_clustering, sep='\t')
paste(commandArgs(), collapse = " ")
print(proc.time() - ptm)
