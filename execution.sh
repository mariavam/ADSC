
#For each dataset is available: rawdata directory, parameters.r, rds objects

script_path='/CEPH/users/mvarea/SingCell_analysis/script'
output_path='/CEPH/users/mvarea/SingCell_analysis/data'
dss=(grubman 
leng
otero
alsema 
mathys)

########################
#DOWNLOAD DATA FROM GEO#
########################
#prefetch ----
#fastq-dump --split-files ----
#cellranger counts ...... ----

########################
#CREATION SEURAT OBJECT#
########################

cd "$output_path/grubman"
Rscript -e "$script_path/prep_grubman.r"

cd "$output_path/leng"
Rscript -e "$script_path/prep_leng.r"

cd "$output_path/otero"
Rscript -e "$script_path/prep_otero.r"

cd "$output_path/alsema"
Rscript -e "$script_path/prep_alsema.r"

cd "$output_path/mathys"
Rscript -e "$script_path/prep_mathys.r"

########################
#QC + FILTERING OF DATA#

#NORMALIZATION#

#PCA#
########################
for DS in dss
do
	cd "$output_path/$DS"
	Rscript -e "$script_path/filter.r"
	Rscript -e "$script_path/normalize.r"
	Rscript -e "$script_path/PCA.r"
done


#############
#INTEGRATION#
#############
cd '$output_path'
Rscript -e "$script_path/integrate.r"


