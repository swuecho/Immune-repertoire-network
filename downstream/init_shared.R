library(rjson) 
library(stringr)

config_file <- "config.json"
config <-fromJSON(file=config_file)

length_cutoff<-as.numeric(config$length_cutoff)
count_cutoff<-as.numeric(config$count_cutoff)
distance_cutoff<-as.numeric(config$distance_cutoff)
patient<-config$patient
reads<-paste('mixcr', config$dbname, sep ='_')

bcr_or_tcr <-config$dbname

compress_name <- function(lang_name) {  
  short_name <-""
  for (i in strsplit(lang_name, "_")) { short_name<-substr(i, 0,1) }
  return(paste(short_name, collapse=""))
}

config_suffix <-paste(config$method, distance_cutoff, 'len', length_cutoff, 'count', count_cutoff, sep='_')

date_str <- format(Sys.time(), "%b_%d_%Y")
#date_str <- "Aug_22_2016"

result_dir <- file.path('./result/', patient, paste(reads,date_str, sep="_"), config_suffix)
# add a method to trim individul part of a string

result_dir <- file.path('./result/', patient, paste(reads,date_str, sep="_"), config_suffix)
                       
config_suffix <- compress_name(config_suffix)

# create target dir
if (!file.exists(result_dir)) {
  dir.create(result_dir, recursive=TRUE)
   # keep a copy of config file
}

file.copy(config_file, result_dir)

raw_csv_file <- file.path('./result/',patient, paste0(patient, '_', bcr_or_tcr ,'.csv'))
raw_data_file <- file.path(result_dir,paste0(patient, '_raw_data_', config_suffix, '.txt'))

data_for_clustering <- file.path(result_dir,paste0(patient, '_cleaned_data_with_majority_vote_', config_suffix, '.txt'))
barcode_count_file <- file.path(result_dir,paste0(patient, '_barcode_and_count_', config_suffix, '.csv'))
subset_count_file <- file.path(result_dir,paste0(patient, '_subset_name_and_cout_', config_suffix, '.csv'))
overall_cluster_edge_list_file <- file.path(result_dir, paste0(patient, "_final_edge_list_of_all_", config_suffix , ".txt"))

all_clusters_csv_file <-file.path(result_dir, paste0(patient, "_cleaned_clusters_all_", config_suffix ,".csv"))
all_clusters_gml_file <- file.path(result_dir, paste0(patient, "_cleaned_clusters_all_", config_suffix, ".gml"))

csf_related_clusters_csv_file <-file.path(result_dir, paste0(patient, "_csf_related_clusters_", config_suffix ,".csv"))
csf_related_clusters_gml_file <- file.path(result_dir, paste0(patient, "_csf_related_clusters_", config_suffix, ".gml"))
