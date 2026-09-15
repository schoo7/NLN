rm(list=ls())

library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
# BiocManager::install("glmGamPoi")
library(glmGamPoi)

seurat_list <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Xiaoling/choushi/diet_seurat_list(240918).rds")

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

s.genes
g2m.genes

#-------------------
celLcycle_pca_pipeline <- function(sobj){
  sobj <- NormalizeData(sobj)
  # Score cells for cell cycle
  sobj <- CellCycleScoring(sobj, g2m.features = g2m.genes, s.features = s.genes)
  sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  sobj <- ScaleData(sobj)
  sobj <- RunPCA(sobj)
  
  return(sobj)
}

# View cell cycle scores and phases assigned to cells                                 
# View(merged_seurat@meta.data)     

for (i in names(seurat_list)){
  sobj <- celLcycle_pca_pipeline(seurat_list[[i]])
  
  DimPlot(sobj,reduction = "pca", group.by= "Phase", split.by = "Phase", raster=FALSE)
  ggsave(paste0("/work/SCCC/s228508/grey_Xiaoling/choushi/cc/cc_phase/",i,".pdf"), width=8, height=6)
  
  DimPlot(sobj, reduction = "pca", group.by= "Phase", raster=FALSE)
  ggsave(paste0("/work/SCCC/s228508/grey_Xiaoling/choushi/cc/cc_phase_2/",i,".pdf"), width=8, height=6)
}






