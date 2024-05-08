######FILE WHERE FROM RAW DATA WE OBTAIN THE SEURAT OBJECT######

#INITIAL DATA AVAILABLE --> Group	Replicate	Region  Sex	Study	Age	PMI

#FINAL METADATA COLUMNS --> orig.ident / Study / Group / Replicate / Sex / Age / PMI / Region

script_path='/CEPH/users/mvarea/SingCell_analysis/scripts/'
output_path='/CEPH/users/mvarea/SingCell_analysis/data/otero/'

source(paste0(script_path,'libraries.r'))
source(paste0(output_path,'parameters.r'))
#source(paste0(script_path,'functions.r'))
###############
#CODE#
###############
if (file.exists(paste0(output_path,"seurat_rawdata.rds"))) {
  stop("The file already exists.")
}

otero.data <- Read10X(data.dir=c(O.AD1A.1 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/1MAP2/outs/filtered_feature_bc_matrix/",
                                 O.AD1B.2 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/1AT8/outs/filtered_feature_bc_matrix/",
                                 O.AD2A.3 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/2MAP2/outs/filtered_feature_bc_matrix/",
                                 O.AD2B.4 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/2AT8/outs/filtered_feature_bc_matrix/",
                                 O.AD3A.5 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/3MAP2/outs/filtered_feature_bc_matrix/",
                                 O.AD3B.6 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/3AT8/outs/filtered_feature_bc_matrix/",
                                 O.AD4A.7 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/4MAP2/outs/filtered_feature_bc_matrix/",
                                 O.AD4B.8 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/4AT8/outs/filtered_feature_bc_matrix/",
                                 O.AD5A.9 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/5MAP2/outs/filtered_feature_bc_matrix/",
                                 O.AD5B.10 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/5AT8/outs/filtered_feature_bc_matrix/",
                                 O.AD6A.11 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/6MAP2/outs/filtered_feature_bc_matrix/",
                                 O.AD6B.12 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/6AT8/outs/filtered_feature_bc_matrix/",
                                 O.AD7A.13 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/7MAP2/outs/filtered_feature_bc_matrix/",
                                 O.AD7B.14 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/7AT8/outs/filtered_feature_bc_matrix/",
                                 O.AD8A.15 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/8MAP2/outs/filtered_feature_bc_matrix/",
                                 O.Ct1A.1 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/CT1/outs/filtered_feature_bc_matrix/",
                                 O.Ct2A.2 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/CT2/outs/filtered_feature_bc_matrix/",
                                 O.Ct3A.3 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/CT3/outs/filtered_feature_bc_matrix/",
                                 O.Ct4A.4 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/CT4/outs/filtered_feature_bc_matrix/",
                                 O.Ct5A.5 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/CT5/outs/filtered_feature_bc_matrix/",
                                 O.Ct6A.6 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/CT6/outs/filtered_feature_bc_matrix/",
                                 O.Ct7A.7 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/CT7/outs/filtered_feature_bc_matrix/",
                                 O.Ct8A.8 = "/CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/CT8/outs/filtered_feature_bc_matrix/"))
                  #O.AD8_B.16 = "CEPH/users/mvarea/SingCell_analysis/data/otero/raw_data/8MAP2/outs/filtered_feature_bc_matrix"

                     
otero.info <- tibble(
  SampleID = c("AD1A","AD1B","AD2A","AD2B","AD3A","AD3B","AD4A","AD4B","AD5A","AD5B","AD6A","AD6B","AD7A","AD7B","AD8A","AD8B","CT1A","CT2A","CT3A","CT4A","CT5A","CT6A","CT7A","CT8A"),
  Age =  c(93,93,79,79,81,81,57,57,89,89,73,73,89,89,62,62,61,67,87,67,72,66,68,71),
  PMI =  c(4,4,19.5,19.5,16,16,14,14,1,1,13,13,1,1,11,11,19.5,11.8,9.3,33,16.4,11.25,19.25,24.9))

############## FILTERING (REMOVE PMI GENES) ############## 
# The nº of PMI genes are 100, because they are the 100 that has past the 5% cutoff (Zhu et al.)
read_csv("/CEPH/users/mvarea/SingCell_analysis/Data/genesPMI.csv",
         col_names = c("pmi.genes"),
         col_types = cols(pmi.genes = col_character())) %>%
  pull(pmi.genes)->genesPMI
# Filtered data: 
otero.data.f <- otero.data[!(rownames(otero.data) %in% genesPMI), ]

############## CREATE SEURAT OBJECT ############## 
# min.cells = 3 --> Sara's and Manual recomend it // min.features = 200 --> Manual recomend it
Otero <- CreateSeuratObject(counts = otero.data.f, project = "Ot_ADSC_prjct", min.cells = goodcells_min, min.features = features_min)

############## EDIT METADATA ############## 
# FINAL METADATA COLUMNS --> orig.ident / Study / Group / Replicate / Sex / Age / PMI / Region

Otero[[]] %>%
  rownames_to_column("cell") %>%
  as_tibble %>%
  separate_wider_delim(orig.ident, delim = ".",
                       names = c("Study","SampleID","Replicate"),
                       cols_remove = FALSE) %>%
  separate_wider_position(SampleID,c(Group = 2, PatientID = 1, Inmuno = 1),
                          cols_remove = FALSE) %>%
  mutate(Study = case_when(Study == "G" ~ "Grubman",
                           Study == "L" ~ "Leng",
                           Study == "O" ~ "Otero",
                           Study == "A" ~ "Alserna",
                           TRUE         ~ NA_character_) %>%
           as.factor,
         Sex = case_when(Group == "AD" ~ case_when(PatientID %in% c(1,2,5,6,7) ~ "F",
                                                   TRUE                     ~ "M"),
                         Group == "Ct" ~ case_when(PatientID %in% c(2,4,7) ~ "F",
                                                   TRUE                  ~ "M")) %>%
           as.factor,
         Region = "PC" %>%
           as.factor,
         Inmuno = case_when(Inmuno == "A" ~ "MAP2+/AT8-",
                            Inmuno == "B" ~ "AT8+") %>%
           as.factor) %>%
  left_join(otero.info, by="SampleID") %>%
  as.data.frame %>%
  column_to_rownames("cell") %>%
  mutate(across(where(is.character), as.factor)) -> Otero

saveRDS(Otero, paste0(output_path,"seurat_rawdata.rds"))
