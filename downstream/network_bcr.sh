#!/bin/bash

patient_id=$(/home/hwu/dotfiles/app/jq .patient config.json)

## to obtain data from mysql database and perform cleaning, QC control and finally write to file
stage="111111111111"
echo $stage
Rscript get_data_bcr.R > /tmp/downstream.log

stage="22222222222"
echo $stage
Rscript edge_list_from_distance.R >> /tmp/downstream.log

stage="333333333"
echo $stage
Rscript write_all_clusters.R >> /tmp/downstream.log
Rscript write_all_clusters_with_subcluster.R >> /tmp/downstream.log


### generate intermediate connctome
stage="444444444"
echo $stage
Rscript connectome_script_all_clsuter_and_write_number_of_clusters_and_write_overlapping_amongsubsets.R  >> /tmp/downstream.log

stage="555555555"
echo $stage
## get simplfied connectome
Rscript connectome_script_CSF_clusters_plus_CSF_to_PB.R >> /tmp/downstream.log

stage="6666666666"
echo $stage
## assess overlap among the subset
Rscript write_subset_overlapping.R >> /tmp/downstream.log

cat /tmp/downstream.log | mailx -s $patient_id echowuhao@gmail.com
cat config.json | mailx -s $patient_id echowuhao@gmail.com

