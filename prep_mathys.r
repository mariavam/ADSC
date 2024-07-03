######FILE WHERE FROM RAW DATA WE OBTAIN THE SEURAT OBJECT######

#INITIAL DATA AVAILABLE --> Group	Replicate	Region  Sex	Study	Age	PMI

#FINAL METADATA COLUMNS --> orig.ident / Study / Group / Replicate / Sex / Age / PMI / Region

script_path='/CEPH/users/mvarea/SingCell_analysis/Code/'
output_path='/CEPH/users/mvarea/SingCell_analysis/data/mathys/'

source(paste0(script_path,'libraries.r'))
source(paste0(script_path,'parameters.r'))
#source(paste0(script_path,'functions.r'))


###############
#CODE#
###############
if (file.exists(paste0(output_path,"prep_mathys.rds"))) {
  stop("The file already exists.")
}

mathys.data <- Read10X(data.dir=c(M.Ct01 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8757/outs/filtered_feature_bc_matrix/"),
                                  M.Ct02 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8761/outs/filtered_feature_bc_matrix/"),
                                  M.Ct03 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8763/outs/filtered_feature_bc_matrix/"),
                                  M.Ct04 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8769/outs/filtered_feature_bc_matrix/"),
                                  M.Ct05 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8773/outs/filtered_feature_bc_matrix/"),
                                  M.Ct06 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8775/outs/filtered_feature_bc_matrix/"),
                                  M.Ct07 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8783/outs/filtered_feature_bc_matrix/"),
                                  M.Ct08 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8785/outs/filtered_feature_bc_matrix/"),
                                  M.Ct09 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8789/outs/filtered_feature_bc_matrix/"),
                                  M.Ct10 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8795/outs/filtered_feature_bc_matrix/"),
                                  M.Ct11 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8777/outs/filtered_feature_bc_matrix/"),
                                  M.Ct12 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8779/outs/filtered_feature_bc_matrix/"),
                                  M.Ct13 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8753/outs/filtered_feature_bc_matrix/"),
                                  M.Ct14 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8755/outs/filtered_feature_bc_matrix/"),
                                  M.Ct15 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8759/outs/filtered_feature_bc_matrix/"),
                                  M.Ct16 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8765/outs/filtered_feature_bc_matrix/"),
                                  M.Ct17 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8767/outs/filtered_feature_bc_matrix/"),
                                  M.Ct18 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8771/outs/filtered_feature_bc_matrix/"),
                                  M.Ct19 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8781/outs/filtered_feature_bc_matrix/"),
                                  M.Ct20 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8787/outs/filtered_feature_bc_matrix/"),
                                  M.Ct21 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8791/outs/filtered_feature_bc_matrix/"),
                                  M.Ct22 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8793/outs/filtered_feature_bc_matrix/"),
                                  M.Ct23 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8797/outs/filtered_feature_bc_matrix/"),
                                  M.Ct24 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8799/outs/filtered_feature_bc_matrix/"),
                                  M.AD01 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8774/outs/filtered_feature_bc_matrix/"),
                                  M.AD02 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8778/outs/filtered_feature_bc_matrix/"),
                                  M.AD03 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8796/outs/filtered_feature_bc_matrix/"),
                                  M.AD04 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8790/outs/filtered_feature_bc_matrix/"),
                                  M.AD05 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8764/outs/filtered_feature_bc_matrix/"),
                                  M.AD06 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8786/outs/filtered_feature_bc_matrix/"),
                                  M.AD07 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8776/outs/filtered_feature_bc_matrix/"),
                                  M.AD08 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8762/outs/filtered_feature_bc_matrix/"),
                                  M.AD09 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8770/outs/filtered_feature_bc_matrix/"),
                                  M.AD10 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8758/outs/filtered_feature_bc_matrix/"),
                                  M.AD11 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8780/outs/filtered_feature_bc_matrix/"),
                                  M.AD12 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8784/outs/filtered_feature_bc_matrix/"),
                                  M.AD13 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8760/outs/filtered_feature_bc_matrix/"),
                                  M.AD14 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8756/outs/filtered_feature_bc_matrix/"),
                                  M.AD15 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8772/outs/filtered_feature_bc_matrix/"),
                                  M.AD16 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8800/outs/filtered_feature_bc_matrix/"),
                                  M.AD17 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8766/outs/filtered_feature_bc_matrix/"),
                                  M.AD18 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8792/outs/filtered_feature_bc_matrix/"),
                                  M.AD19 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8768/outs/filtered_feature_bc_matrix/"),
                                  M.AD20 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8782/outs/filtered_feature_bc_matrix/"),
                                  M.AD21 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8754/outs/filtered_feature_bc_matrix/"),
                                  M.AD22 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8794/outs/filtered_feature_bc_matrix/"),
                                  M.AD23 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8798/outs/filtered_feature_bc_matrix/"),
                                  M.AD24 = c("/CEPH/users/mvarea/SingCell_analysis/data/mathys/raw_data/D17-8788/outs/filtered_feature_bc_matrix/")))

mathys.info <- tibble(
  SampleID=c("Ct01","Ct02","Ct03","Ct04","Ct05","Ct06","Ct07","Ct08","Ct09","Ct10","Ct11","Ct12","Ct13","Ct14","Ct15","Ct16","Ct17","Ct18","Ct19","Ct20","Ct21","Ct22","Ct23","Ct24",
         "AD01","AD02","AD03","AD04","AD05","AD06","AD07","AD08","AD09","AD10","AD11","AD12","AD13","AD14","AD15","AD16","AD17","AD18","AD19","AD20","AD21","AD22","AD23","AD24"),
  PMI=c(4.16,15.3,18.2,4.5,7,5.3,3.25,3.6,16.72,6.3,3.9,20.8,
        1.33,8.6,16.8,15.17,3.25,2.6,
        16.7,23,6.5,7.6,85,2.3,
        3.3,4.3,18.2,26.2,14,NA,
        2.2,4.3,7.8,5.1,16.7,17,
        20.5,18,13.3,1.5,4.5,4.5,
        4.3,4.6,13.8,8.7,2.3,1.5),
  Age=c(90,90,79,90,89,84,90,87,87.3,90,90,82,80,88,78.2,85,77,90,81,80,86,83,87,87,
        81,90,90,88,76,88,83,90,90,75,84,90,86.4,88.5,89,85.8,85.8,87,84,80.5,89,83.5,87,86.3)
)

############## FILTERING (REMOVE PMI GENES) ############## 
# The nº of PMI genes are 100, because they are the 100 that has past the 5% cutoff (Zhu et al.)
read_csv("/CEPH/users/mvarea/SingCell_analysis/data/genesPMI.csv",
         col_names = c("pmi.genes"),
         col_types = cols(pmi.genes = col_character())) %>%
  pull(pmi.genes)->genesPMI
# Filtered data: 
mathys.data.f <- mathys.data[!(rownames(mathys.data) %in% genesPMI), ]

############## CREATE SEURAT OBJECT ############## 
# min.cells = 3 --> Sara's and Manual recomend it // min.features = 200 --> Manual recomend it
Mathys <- CreateSeuratObject(counts = mathys.data.f, project = "Mat_ADSC_prjct", min.cells = goodcells_min, min.features = features_min)

############## EDIT METADATA ############## 
# FINAL METADATA COLUMNS --> orig.ident / Study / Group / Replicate / Sex / Age / PMI / Region

Mathys[[]] %<>%
  rownames_to_column("cell") %>%
  as_tibble %>%
  separate_wider_delim(orig.ident, delim = ".",
                       names = c("Study","SampleID"),
                       cols_remove = FALSE) %>%
  separate_wider_position(SampleID,c(Group = 2, Replicate = 2),
                          cols_remove = FALSE) %>%
  mutate(Study = case_when(Study == "G" ~ "Grubman",
                           Study == "L" ~ "Leng",
                           Study == "O" ~ "Otero",
                           Study == "A" ~ "Alserna",
                           Study == "M" ~ "Mathys",
                           TRUE         ~ NA_character_) %>%
           as.factor,
         Sex = case_when(Replicate>12 ~ "M",
                         TRUE         ~ "F") %>%
           as.factor,
         Region = "DPC" %>%
           as.factor) %>%
  left_join(mathys.info, by="SampleID") %>%
  as.data.frame %>%
  column_to_rownames("cell") %>%
  mutate(across(where(is.character), as.factor))

saveRDS(Mathys, paste0(output_path,"seurat_rawdata.rds"))
