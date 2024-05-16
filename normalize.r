script_path='/CEPH/users/mvarea/SingCell_analysis/scripts/'
output_path='/CEPH/users/mvarea/SingCell_analysis/data'

source(paste0(script_path,'libraries.r'))
source('parameters.r')
#source(paste0(script_path,'functions.r'))



if (file.exists("seurat_normalized.rds")) {
  stop("The file already exists.")
}

########### NORMALIZATION
Data <- readRDS("seurat_filtered.rds")

Data[["RNA"]] <- split(Data[["RNA"]], f=Data$orig.ident)

#Now we perform normally the SCTransform function
Data <- SCTransform(Data, variable.features.n =nfeatures, ncells = ncells,
                    vars.to.regress = c("nFeature_RNA", "percent.mt"), 
                    seed.use = seeduse, verbose = FALSE, vst.flavor = "v2")

saveRDS(Data,"seurat_normalized.rds")