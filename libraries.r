library(Seurat) #Toolkit for quality control, analysis and exploration of scRNA seq data. 
library(gdata) #Data manipulation (combining objects) --> used by Sara bfore create the Study column into meta.data. 
library(readr) #Reed quickly CSV/TSV easy use 
library(dplyr) #Manipulate data intuitively and efficiently Filter/Select/Modify/Resume %>%
library(tidyr) #Order and manipulate data -> variable = column / observ = rows %>%
library(tibble) # For rownames_to_column etc
library(ggplot2) #Generate graphs (PDF guide)
library(plotly) #Generate graphs
library(magrittr) #Chaining commands with %>%
library(clustree)
library(future)
library(cowplot)
library(patchwork)