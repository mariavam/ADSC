######FILE WHERE FROM RAW DATA WE OBTAIN THE SEURAT OBJECT######

#INITIAL DATA AVAILABLE --> Group	Replicate	Region  Sex	Study	Age	PMI

#FINAL METADATA COLUMNS --> orig.ident / Study / Group / Replicate / Sex / Age / PMI / Region

script_path='/CEPH/users/mvarea/SingCell_analysis/Code/'
output_path='/CEPH/users/mvarea/SingCell_analysis/data/grubman/'

source(paste0(script_path,'libraries.r'))
source(paste0(output_path,'parameters.r'))
#source(paste0(script_path,'functions.r'))

###############
#CODE#
###############
if (file.exists(paste0(output_path,"seurat_rawdata.rds"))) {
  stop("The file already exists.")
}

grub.data <- Read10X(data.dir=c(G.AD1="/CEPH/users/mvarea/SingCell_analysis/data/grubman/raw_data/A1A2/",
                                G.AD2="/CEPH/users/mvarea/SingCell_analysis/data/grubman/raw_data/A3A4/",
                                G.AD3="/CEPH/users/mvarea/SingCell_analysis/data/grubman/raw_data/A5A6/",
                                G.Ct1="/CEPH/users/mvarea/SingCell_analysis/data/grubman/raw_data/C1C2/",
                                G.Ct2="/CEPH/users/mvarea/SingCell_analysis/data/grubman/raw_data/C3C4/",
                                G.Ct3="/CEPH/users/mvarea/SingCell_analysis/data/grubman/raw_data/C5C6/"))


Grub.info <- tibble(
  Case = c("AD1.1", "AD1.2","AD2.3","AD2.4","AD3.5","AD3.6","Ct1.1","Ct1.2","Ct2.3","Ct2.4","Ct3.5","Ct3.6"),
  Age =  c(91,83.8,67.8,83,73,74.6,67.3,82.7,72.6,75.6,77.5,82.7),
  PMI =  c(8,10,21,34,9.5,30,24,28.5,42.5,46,53.5,27))

############## FILTERING (REMOVE PMI GENES) ############## 
# The nº of PMI genes are 100, because they are the 100 that has past the 5% cutoff (Zhu et al.)
read_csv("/CEPH/users/mvarea/SingCell_analysis/data/genesPMI.csv",
         col_names = c("pmi.genes"),
         col_types = cols(pmi.genes = col_character())) %>%
  pull(pmi.genes)->genesPMI
# Filtered data: 
grub.data.f <- grub.data[!(rownames(grub.data) %in% genesPMI), ]

############## CREATE SEURAT OBJECT ############## 
# min.cells = 3 --> Sara's and Manual recomend it // min.features = 200 --> Manual recomend it
Grub <- CreateSeuratObject(counts = grub.data.f, project = "Grub_ADSC_prjct", min.cells = goodcells_min, min.features = features_min)

############## EDIT METADATA ############## 
# FINAL METADATA COLUMNS --> orig.ident / Study / Group / Replicate / Sex / Age / PMI / Region

Grub[[]] %<>%
  rownames_to_column("cell") %>%
  as_tibble %>%
  separate_wider_delim(orig.ident, delim = ".",
                       names = c("Study","Sample"),
                       cols_remove = FALSE) %>%
  separate_wider_position(Sample, c(Group = 2, Replicate = 1)) %>%
  mutate(Study = case_when(Study == "G" ~ "Grubman",
                           Study == "L" ~ "Leng",
                           Study == "O" ~ "Otero",
                           Study == "A" ~ "Alserna",
                           TRUE         ~ NA_character_) %>%
           as.factor,
         Sex = case_when(orig.ident %in% c("G.AD2", "G.Ct1") ~ "F",
                         TRUE                                ~ "M") %>%
           as.factor,
         Region = "EC" %>%
           as.factor) %>%
  left_join(Grub.info %>%
              separate_wider_delim(Case,".",names = c("Sample","Replicate")) %>%
              group_by(Sample) %>%
              summarise(across(Age:PMI, mean), .groups = "drop") %>%
              (dplyr::rename)(orig.ident=Sample) %>%
              mutate(orig.ident = paste0("G.", orig.ident))) %>%
  as.data.frame %>%
  column_to_rownames("cell") %>%
  mutate(across(where(is.character), as.factor))

saveRDS(Grub, paste0(output_path,"seurat_rawdata.rds"))
