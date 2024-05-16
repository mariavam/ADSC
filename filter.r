
script_path='/CEPH/users/mvarea/SingCell_analysis/scripts/'
output_path='/CEPH/users/mvarea/SingCell_analysis/data'

source(paste0(script_path,'libraries.r'))
source('parameters.r')
#source(paste0(script_path,'functions.r'))



if (file.exists("seurat_filtered.rds")) {
  stop("The file already exists.")
}

#############QC 
Data <- readRDS("seurat_rawdata.rds")

Data[["percent.mt"]] <- PercentageFeatureSet(Data, pattern = "^MT-")
saveRDS(Data, "seurat_graphsfiltering.rds")

##############Filtering 1
#Filtering 2: Limits determinated by Visualization (07 & 08)
Data <- subset(Data, subset = nFeature_RNA > nFeat_min & nFeature_RNA < nFeat_max & percent.mt < perc_limit)

'
#Statistics
discarted_cells <- nrow(Data[[]]) - nrow(Data2[[]])
percentage_discells <- (discarted_cells / nrow(Data[[]]))*100
percentage_filtcells <- (nrow(Data2) / nrow(Data[[]]))*100
'
#Filtering 2: Using minimum number good cells. 
m <- GetAssayData(Data, layer = "counts")
t <- table(m@i)
head(t, 20)

t.selected <- t[t >= goodcells_min]
head(t.selected, 20)

genes.selected <- rownames(m)[as.numeric(names(t.selected)) + 1]
length(genes.selected)

Data<- subset(Data, features = genes.selected)
Data

saveRDS(Data,"seurat_filtered.rds")
        
        
        
