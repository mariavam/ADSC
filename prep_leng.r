######FILE WHERE FROM RAW DATA WE OBTAIN THE SEURAT OBJECT######

#INITIAL DATA AVAILABLE -->

#FINAL METADATA COLUMNS --> orig.ident / Study / Group / Replicate / Sex / Age / PMI / Region
script_path='/CEPH/users/mvarea/SingCell_analysis/Code/'
output_path='/CEPH/users/mvarea/SingCell_analysis/data/leng/'

source(paste0(script_path,'libraries.r'))
source(paste0(output_path,'parameters.r'))
#source(paste0(script_path,'functions.r'))
###############
#CODE#
###############
if (file.exists(paste0(output_path,"seurat_rawdata.rds"))) {
  stop("The file already exists.")
}

leng.data <- Read10X(data.dir=c(L.SFG.AD1 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/SFG1/outs/filtered_feature_bc_matrix/.'),
                                L.SFG.AD2 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/SFG2/outs/filtered_feature_bc_matrix/.'),
                                L.SFG.AD3 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/SFG3/outs/filtered_feature_bc_matrix/.'),
                                L.SFG.AD4 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/SFG4/outs/filtered_feature_bc_matrix/.'),
                                L.SFG.AD5 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/SFG5/outs/filtered_feature_bc_matrix/.'),
                                L.SFG.Ct1 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/SFG6/outs/filtered_feature_bc_matrix/.'),
                                L.SFG.Ct2 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/SFG7/outs/filtered_feature_bc_matrix/.'),
                                L.SFG.Ct3 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/SFG8/outs/filtered_feature_bc_matrix/.'),
                                L.SFG.Ct4 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/SFG9/outs/filtered_feature_bc_matrix/.'),
                                L.SFG.Ct5 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/SFG10/outs/filtered_feature_bc_matrix/.'),
                                L.EC.AD1 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/EC1/outs/filtered_feature_bc_matrix/.'),
                                L.EC.AD2 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/EC2/outs/filtered_feature_bc_matrix/.'),
                                L.EC.AD3 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/EC3/outs/filtered_feature_bc_matrix/.'),
                                L.EC.AD4 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/EC4/outs/filtered_feature_bc_matrix/.'),
                                L.EC.AD5 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/EC5/outs/filtered_feature_bc_matrix/.'),
                                L.EC.Ct1 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/EC6/outs/filtered_feature_bc_matrix/.'),
                                L.EC.Ct2 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/EC7/outs/filtered_feature_bc_matrix/.'),
                                L.EC.Ct3 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/EC8/outs/filtered_feature_bc_matrix/.'),
                                L.EC.Ct4 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/EC9/outs/filtered_feature_bc_matrix/.'),
                                L.EC.Ct5 = c('/CEPH/users/mvarea/SingCell_analysis/data/leng/raw_data/EC10/outs/filtered_feature_bc_matrix/.')))

leng.info <- tibble(Agevalues = c(50, 60, 71, 72, 77, 87, 91, 72, 82, 82),
                    PMIvalues = c(13, 12, 12, 15, 4.9, 30, 50, 6.9, 6.7, 9))

############## FILTERING (REMOVE PMI GENES) ############## 
# The nº of PMI genes are 100, because they are the 100 that has past the 5% cutoff (Zhu et al.)
read_csv("/CEPH/users/mvarea/SingCell_analysis/data/genesPMI.csv",
         col_names = c("pmi.genes"),
         col_types = cols(pmi.genes = col_character())) %>%
         pull(pmi.genes)->genesPMI
# Filtered data: 
leng.data.f <- leng.data[!(rownames(leng.data) %in% genesPMI), ]

############## CREATE SEURAT OBJECT ############## 
# min.cells = 3 --> Saras and Manual recomend it // min.features = 200 --> Manual recomend it
Leng <- CreateSeuratObject(counts = leng.data.f, project = "Leng_ADSC_prjct", min.cells = goodcells_min, min.features = features_min)

############## EDIT METADATA ############## 
# FINAL METADATA COLUMNS --> orig.ident / Study / Group / Replicate / Sex / Age / PMI / Region
Agevalues = c(50, 60, 71, 72, 77, 87, 91, 72, 82, 82)
PMIvalues = c(13, 12, 12, 15, 4.9, 30, 50, 6.9, 6.7, 9)

Leng[[]] %<>%
  rownames_to_column("cell") %>%
  as_tibble %>%
  separate_wider_delim(orig.ident, delim = ".",
                       names = c("Study","Region","SampleID"),
                       cols_remove = FALSE) %>%
  separate_wider_position(SampleID,c(Group = 2, PatientID = 1),
                          cols_remove = FALSE) %>%
  mutate(Study = case_when(Study == "G" ~ "Grubman",
                           Study == "L" ~ "Leng",
                           Study == "O" ~ "Otero",
                           Study == "A" ~ "Alserna",
                           TRUE         ~ NA_character_) %>%
           as.factor,
         Sex = "M" %>%
           as.factor,
         Age = ifelse(PatientID %in% 1:10, Agevalues[as.integer(PatientID)], NA) %>%
           as.factor,
         PMI = ifelse(PatientID %in% 1:10, PMIvalues[as.integer(PatientID)], NA) %>%
           as.factor) %>%
  left_join(arrange(., PatientID) %>%
              group_by(Group) %>%
              distinct(PatientID) %>% #Makesure in Replicate 1 Replicate <-> 1 Patient ID
              mutate(Replicate =as.factor(1:n())) %>% #Numerate the samples in each group (Ct 1 2 3 / AD 1 2 3 4...)
              ungroup) %>%  as.data.frame %>%
  column_to_rownames("cell") %>%
  mutate(across(where(is.character), as.factor))

saveRDS(Leng, paste0(output_path,"seurat_rawdata.rds"))


' OTHER OPTION
SFG_object <- readRDS("/CEPH/users/mvarea/SingCell_analysis/Data/Leng_data/sce.SFG.scAlign.assigned.rds") %>%
  as.Seurat
EC_object <- readRDS("/CEPH/users/mvarea/SingCell_analysis/Data/Leng_data/sce.EC.scAlign.assigned.rds") %>%
  as.Seurat

Leng <- merge(SFG_object, EC_object)

Agevalues = c(50, 60, 71, 72, 77, 87, 91, 72, 82, 82)
PMIvalues = c(13, 12, 12, 15, 4.9, 30, 50, 6.9, 6.7, 9)

Leng[[]]%>%
  select(Braak = BraakStage,
         PatientID, 
         ClusterLeng = clusterAssignment,
         Region = BrainRegion, #Select columns + Change names to both columns: Braak and ClusterLeng
  ) %>%
  mutate(Group = ifelse(Braak %in% 0, "Ct", "AD") %>% #Create column Group
           factor(levels = c("Ct", "AD"))) %>% #Create columns orig.ident and Sex) 
  left_join(arrange(., PatientID) %>%
              group_by(Group) %>%
              distinct(PatientID) %>% #Makesure in Replicate 1 Replicate <-> 1 Patient ID
              mutate(Replicate =as.factor(1:n())) %>% #Numerate the samples in each group (Ct 1 2 3 / AD 1 2 3 4...)
              ungroup) %>%
  mutate(Age = ifelse(PatientID %in% 1:10, Agevalues[as.integer(PatientID)], NA),
         PMI = ifelse(PatientID %in% 1:10, PMIvalues[as.integer(PatientID)], NA),
         orig.ident = paste0("L.",Region,".",Group, Replicate),
         Sex = factor("M", levels = c("F", "M")) %>%
           as.factor) %>%
  relocate(orig.ident, Group, Replicate) %>% 
  CreateSeuratObject(counts = GetAssayData(Leng, "originalexp", "counts"), #Indicate to create a new Seurat Object from one of the layers of the first one.
                     meta.data = .) -> Leng

Leng[["Braak"]] <- NULL
Leng[["PatientID"]] <- NULL
Leng[["ClusterLeng"]] <- NULL
'