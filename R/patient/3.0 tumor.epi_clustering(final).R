rm(list=ls())
library(RColorBrewer)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)

seurat_list <- readRDS("/work/SCCC/s228508/choushi/epi_seurat_list.rds")
tumor.epi <- merge(x = seurat_list[[1]], y = seurat_list[-1])
gene_list <- readRDS("/work/SCCC/s228508/choushi/gene_list_wnt.rds")

#-------------------------------------------------------------------------
# wnt_features <- c()
# for (i in (1:4)) {
#   wnt_features <- append(wnt_features, gene_list[[i]])
# }
# 
# wnt_features <- unique(wnt_features)
# 
# DefaultAssay(tumor.epi) <- "RNA"
# tumor.epi <- NormalizeData(tumor.epi)
# tumor.epi <- ScaleData(tumor.epi,features = wnt_features, verbose = T)
# tumor.epi <- RunPCA(tumor.epi, features = wnt_features, npcs = 50, verbose = T)

#-------------------------------------------------------------------------
features_list <- read.table("/work/SCCC/s228508/choushi/top2000.csv")
# features <- gene_list$genes
features <- head(features_list$genes , 1000)

DefaultAssay(tumor.epi) <- "RNA"
tumor.epi <- NormalizeData(tumor.epi)
# tumor.epi <- FindVariableFeatures(tumor.epi, selection.method = "vst",nfeatures = 2000)
tumor.epi <- ScaleData(tumor.epi,features = features, verbose = T)
tumor.epi <- RunPCA(tumor.epi, features = features, npcs = 50, verbose = T)

#-------------------------------------------------------------------------
# DefaultAssay(tumor.epi) <- "RNA" 
# tumor.epi <- NormalizeData(tumor.epi)
# tumor.epi <- FindVariableFeatures(tumor.epi, selection.method = "vst",nfeatures = 2000)
# tumor.epi <- ScaleData(tumor.epi, verbose = T) 
# tumor.epi <- RunPCA(tumor.epi, npcs = 50, verbose = T)

#-----------------------------
DimPlot(object = tumor.epi, reduction = "pca", group.by = "Patient_ID")
# VlnPlot(object = tumor.epi, features = "PC_1", group.by = "Patient_ID")

#-----------------------------
### integration
tumor.epi <- RunHarmony(tumor.epi, group.by.vars = "Patient_ID",lambda = 1)

DimPlot(object = tumor.epi, reduction = "harmony", group.by = "Patient_ID")
# VlnPlot(object = tumor.CD8_harmony, features = "harmony_1", group.by = "Patient_ID")
ElbowPlot(tumor.epi, ndims = 50)

#--------------------------------------------------------------------------------------------
set.seed(123)
pcSelect = 20
tumor.epi <- FindNeighbors(tumor.epi, reduction = "harmony", dims = 1:pcSelect)
tumor.epi <- FindClusters(tumor.epi, resolution = 0.1) 
tumor.epi <- RunUMAP(tumor.epi, reduction = "harmony", dims = 1:pcSelect)

# saveRDS(tumor.epi, "/work/SCCC/s228508/choushi/tumor.epi_top1000_p20_res0.1.rds")

DimPlot(tumor.epi, label = T, raster= F) 

DimPlot(tumor.epi, label = T, raster = F, split.by = "Patient_ID")  

FeaturePlot(object = tumor.epi, features = c("NLN"), pt.size = 1,order = T, label = F,cols = c("#FFF1ED", "#3F0A3F"))
FeaturePlot(object = tumor.epi, features = "NLN", pt.size = 1, order = T, label = T)
VlnPlot(object = tumor.epi, features = "NLN")
DotPlot(object = tumor.epi, features = "NLN", group.by = "Patient_ID")
VlnPlot(object = tumor.epi, features = "NLN", group.by = "Patient_ID")


FeaturePlot(object = tumor.epi, features = "NLN",cols = c("#FFF1ED", "#3F0A3F"))



#---------------------------------------
for(gn in names(gene_list)){
  tumor.epi <- AddModuleScore(tumor.epi, 
                              features = list(gene_list[[gn]]), 
                              assay = "RNA",
                              name = gn)
}

DotPlot(object = tumor.epi, features = paste(names(gene_list), "1", sep=""))+ RotatedAxis()
FeaturePlot(object = tumor.epi, features = paste(names(gene_list), "1", sep=""),cols = c("white", "#3F0A3F"))+ RotatedAxis()




#-------------------------------------------------------------------------------------------------------------------
# save figs
for(gn in names(gene_list))
{
  FeaturePlot(object=tumor.epi, features=paste(gn, "1", sep=""), cols = c("grey", "red"),
              label = TRUE)
  ggsave(paste0("/work/SCCC/s228508/subclustering_figs/feature_plots_pc30_res0.3/",
                gn,".pdf"), width=5, height = 5)
}


showGenes = c("CHGA","CHGB","ENO2","SYP","ENO1","ENO2","ASCL1","PCSK1","CHRNB2","SCG3","SCN3A","ELAVL4","NKX2-1") 
pdf(file="/work/SCCC/s228508/reclustering/draft/showgenes.pdf",width=13,height=8)
FeaturePlot(object = tumor.epi, features = showGenes, cols = c("grey", "red"))
dev.off()

