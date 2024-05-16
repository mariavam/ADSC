###### ! common in all datasets
goodcells_min=3 # --> Criteria followed in the manual (probably statistical robustness / biological relevance / Practical considerations / Consistency)
features_min=200 # --> Manual recomend it
#####specific in each dataset#####
#Grubman
G_nFeat_max=6000
G_nFeat_min=250
G_perc_limit=5

#Leng
L_nFeat_max=9000
L_nFeat_min=300
L_perc_limit=10

#Otero
O_nFeat_max=7000
O_nFeat_min=400
O_perc_limit=10
  
#Alsema
A_nFeat_max=2000
A_nFeat_min=300
A_perc_limit=5

###### !
nfeatures = 3000
ncells = 5000
seeduse = 12345


VDL_dim = 5

heat_cells = 500
heat_dims = 2

elbow_ndims = 40

umap_dims = 30

resol = 0.5
num_cluster = 2
filter_avglog2FC = 1
threshold_logfc = 0.25
slice_head = 10
