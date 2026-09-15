rm(list=ls())
gc()

library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

rds_files_names <- dir("/work/SCCC/s228508/grey_Mu_lab/choushi/after_qc/")
rds_files_names
dir <- paste0("/work/SCCC/s228508/grey_Mu_lab/choushi/after_qc/", rds_files_names)

samples_name <- c("sgNLN_KIF","sgNLN_WNT", "sgNLN", "sgNT")

seurat_list <- list()
for (i in 1:length(dir)) {
  seurat_list[[i]] <- readRDS(dir[i])
  # add refix for each cell barcode to avoid repeats after merge
  seurat_list[[i]] <- RenameCells(seurat_list[[i]], add.cell.id = samples_name[i])
}

names(seurat_list) <- samples_name

#----------------------------------------------------------------------------------------------------------------
for (seurat_obj_name in names(seurat_list)) {
  seurat_list[[seurat_obj_name]]$Sample_ID <- seurat_obj_name
}


# seurat_list[["RPM_173_284M"]]$condition <- "regional_met"
# seurat_list[["RPMAN_290_302M"]]$condition <- "regional_met"
# seurat_list[["RPMA_MALM"]]$condition <- "regional_met"
# seurat_list[["RPMN_MNLM"]]$condition <- "regional_met"

saveRDS(seurat_list, "/work/SCCC/s228508/grey_Mu_lab/choushi/seurat_list(240918).rds")
# 
# #---------------------------------------------------------------------------------------------------------------
# write.table(a,
#           file = "/work/SCCC/s228508/integration/integration_genes_3004.csv",
#           row.names = F,
#           col.names = "integration_gene")


#########################################################################################################################
#########################################################################################################################
seurat_list <- readRDS("/work/SCCC/s228508/grey_Mu_lab/choushi/seurat_list(240918).rds") 

# diet seurat objects
for (seurat_obj_name in names(seurat_list)) {
  seurat_list[[seurat_obj_name]]@commands <- list() 
  
  diet_obj <- CreateSeuratObject(counts = seurat_list[[seurat_obj_name]]@assays[["RNA"]]@counts)
  # Copy relevant metadata columns
  metadata_subset <- seurat_list[[seurat_obj_name]]@meta.data[, c("nCount_RNA", "nFeature_RNA","percent.mt","Sample_ID")]
  rownames(metadata_subset) <- rownames(seurat_list[[seurat_obj_name]]@meta.data)
  diet_obj@meta.data <- metadata_subset
  
  saveRDS(diet_obj, paste0("/work/SCCC/s228508/grey_Mu_lab/choushi/diet_seurat/",seurat_obj_name,".rds"))
}

#---------------------------------------------------------------------------------------------------------------------------
rds_files_names <- dir("/work/SCCC/s228508/grey_Mu_lab/choushi/diet_seurat2/")
rds_files_names
dir <- paste0("/work/SCCC/s228508/grey_Mu_lab/choushi/diet_seurat2/", rds_files_names)

samples_name <- c("sgNLN_WNT", "sgNLN", "sgNT")

diet_seurat_list <- list()
for (i in 1:length(dir)) {
  diet_seurat_list[[i]] <- readRDS(dir[i])
}

names(diet_seurat_list) <- samples_name

saveRDS(diet_seurat_list, "/work/SCCC/s228508/grey_Mu_lab/choushi/diet_seurat_list(241021).rds")











