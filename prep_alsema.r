######FILE WHERE FROM RAW DATA WE OBTAIN THE SEURAT OBJECT######

#INITIAL DATA AVAILABLE --> Group	Replicate	Region  Sex	Study	Age	PMI

#FINAL METADATA COLUMNS --> orig.ident / Study / Group / Replicate / Sex / Age / PMI / Region

script_path='/CEPH/users/mvarea/SingCell_analysis/Code/'
output_path='/CEPH/users/mvarea/SingCell_analysis/Output/'

source(paste0(script_path,'libraries.r'))
source(paste0(script_path,'parameters.r'))
#source(paste0(script_path,'functions.r'))


###############
#CODE#
###############
if (file.exists(paste0(output_path,"seurat_rawdata.rds"))) {
  stop("The file already exists.")
}

alsema.data <- Read10X(data.dir=c(A.AD1 = "/CEPH/users/mvarea/SingCell_analysis/data/alsema/AD1/outs/filtered_feature_bc_matrix/",
                                  A.AD2 = "/CEPH/users/mvarea/SingCell_analysis/data/alsema/AD2/outs/filtered_feature_bc_matrix/"))
alsema.info <- tibble(
  SampleID = c("AD1","AD2"),
  Age = c(76,81),
  PMI = c(3.45,4.05),
  PatientID = c("2019-010","2018-135"),
  Replicate = c(1,2))

############## FILTERING (REMOVE PMI GENES) ############## 
# The nº of PMI genes are 100, because they are the 100 that has past the 5% cutoff (Zhu et al.)
read_csv("/CEPH/users/mvarea/SingCell_analysis/Data/genesPMI.csv",
         col_names = c("pmi.genes"),
         col_types = cols(pmi.genes = col_character())) %>%
  pull(pmi.genes)->genesPMI
# Filtered data: 
alsema.data.f <- alsema.data[!(rownames(alsema.data) %in% genesPMI), ]

############## CREATE SEURAT OBJECT ############## 
# min.cells = 3 --> Sara's and Manual recomend it // min.features = 200 --> Manual recomend it
Alsema <- CreateSeuratObject(counts = alsema.data.f, project = "Als_ADSC_prjct", min.cells = goodcells_min, min.features = features_min)

############## EDIT METADATA ############## 
# FINAL METADATA COLUMNS --> orig.ident / Study / Group / Replicate / Sex / Age / PMI / Region

Alsema[[]] %>%
  rownames_to_column("cell") %>%
  as_tibble %>%
  separate_wider_delim(orig.ident, delim = ".",
                       names = c("Study","SampleID"),
                       cols_remove = FALSE) %>%
  separate_wider_position(SampleID,c(Group = 2, PatientID = 1),
                          cols_remove = FALSE) %>%
  mutate(Study = case_when(Study == "G" ~ "Grubman",
                           Study == "L" ~ "Leng",
                           Study == "O" ~ "Otero",
                           Study == "A" ~ "Alserna",
                           TRUE         ~ NA_character_) %>%
           as.factor,
         Sex = "F" %>%
           as.factor,
         Region = "LPS" %>%
           as.factor) %>%
  left_join(alsema.info, by="SampleID") %>%
  as.data.frame %>%
  column_to_rownames("cell") %>%
  mutate(across(where(is.character), as.factor))->Alsema

saveRDS(Alsema, paste0(output_path,"seurat_rawdata.rds"))
