#!/bin/bash

script_path='/CEPH/users/mvarea/SingCell_analysis/Code'
output_path='/CEPH/users/mvarea/SingCell_analysis/Output/Leng'


#PARAMETERS  --> $FILE_PATH/parameters.r

#LIBRARIES --> $FILE_PATH/libraries.r

#FUNCTIONS --> $FILE_PATH/functions.r


########################
#DOWNLOAD DATA FROM GEO#
########################

# Download Raw Data from GEO platform. The SRA numbers are save with the
# metadata information of each study in Data_information.csv

chmod +x $script_path/download.sh #change permissions for execution

$script_path/download.sh


##################################
#CREATION INITIAL SEURAT OBJECTS#
##################################

# In this step, I'm going to create a seurat object from each dataset
# separatedly. Here Im going to edit them with the properly data. 

# Common metadata information: orig.ident / study / group / replicate / sex /
# age / region /
# The metadata of each study is saved in Data_information.csv, as well as
# initial comments in each R script

Rscript -e "$script_path/prep_grubman.r"

Rscript -e "$script_path/prep_leng.r"

Rscript -e "$script_path/prep_otero.r"

Rscript -e "$script_path/prep_alsema.r"

Rscript -e "$script_path/prep_mathys.r"


#####################
#PREPROCESSING DATA#
#####################


###############################
#PROCESSING ALL DATA TOGUETHER#
###############################
