rm(list=ls())

library(Seurat)
library(flexmix)
library(ggplot2)
library(patchwork)
library(BPCells)
library(SeuratObject)
# library(SeuratData)
library(SeuratWrappers)
library(sctransform)
library(DoubletFinder)
library(fishpond)
library(SummarizedExperiment)
library(UCell)
library(GSEABase)
library(escape)
library(SCpubr)
library(decoupleR)
library(ggsankey)
library(extrafont)
library(extrafontdb)
library(SeuratExtend)

options(Seurat.object.assay.version = "v5")

hupsa=readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/HuPSA_share.rds")
DimPlot2(hupsa, raster = F, theme = theme_umap_arrows())

table(hupsa$sample)
table(hupsa$histo)
length(unique(hupsa$sample))

#-----------------------
Idents(hupsa) <- "histo"
to_remove <- c("Benign", "Normal", "Normal_adj","PCa_Cribriform")  
hupsa <- subset(hupsa, idents = to_remove, invert = TRUE)
# DimPlot(hupsa, raster = F)

Idents(hupsa) <- "cell_type"
to_remove <- c("Normal","Basal","Club","KRT7","Fibroblast","SM","Pericyte","Endothelial","T_CD4","T_CD8","T_proliferating","NK","B",
               "Plasma","Mast","Monocyte","Macrophage","Dendritic")  
hupsa <- subset(hupsa, idents = to_remove, invert = TRUE)
# DimPlot(hupsa, raster = F)
length(unique(hupsa$sample))

##########################
DefaultAssay(hupsa) <- "RNA"
hupsa <- NormalizeData(hupsa)
hupsa <- FindVariableFeatures(hupsa, nfeatures = 200) #200
hupsa <- ScaleData(hupsa)  
hupsa <- RunPCA(hupsa, npcs = 50, verbose = TRUE)
# hupsa <- RunHarmony(hupsa, group.by.vars = "sample")
ElbowPlot(hupsa, ndims = 50)

set.seed(123)
pcSelect <- 20 #10
hupsa <- FindNeighbors(hupsa, reduction = "pca", dims = 1:pcSelect)
hupsa <- FindClusters(hupsa, resolution = 0.08, cluster.name = "clusters_0.1")
hupsa <- RunUMAP(hupsa,reduction = "pca",dims = 1:pcSelect, reduction.name = "umap",return.model = T)

# saveRDS(hupsa, "/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/HuPSA_tumoronly(250430).rds")

DimPlot2(hupsa, group.by = "seurat_clusters", label = T,raster = F,theme = theme_umap_arrows())
DimPlot2(hupsa, group.by = "cell_type", label = T,raster = F,theme = theme_umap_arrows())
DimPlot2(hupsa, split.by = "histo", label = T,raster = F,theme = theme_umap_arrows())

SCpubr::do_DimPlot(sample = hupsa,pt.size = 0.1,raster = F,border.size = 20,legend.position = "none",label = T,label.color = "black")
### data clean


###########################################
# 2nd
# remove some clusters
Idents(hupsa) <- "seurat_clusters"
to_remove <- c("10","11")  
hupsa2 <- subset(hupsa, idents = to_remove, invert = TRUE)
DimPlot(hupsa2, raster = F)

# DefaultAssay(hupsa2) <- "RNA"
# hupsa2 <- NormalizeData(hupsa2)
# hupsa2 <- FindVariableFeatures(hupsa2, nfeatures = 100)
# hupsa2 <- ScaleData(hupsa2)  
# hupsa2 <- RunPCA(hupsa2, npcs = 50, verbose = TRUE)
# # hupsa <- RunHarmony(hupsa, group.by.vars = "sample")
# ElbowPlot(hupsa2, ndims = 50)

set.seed(123)
pcSelect <- 25
hupsa2 <- FindNeighbors(hupsa2, reduction = "pca", dims = 1:pcSelect)
hupsa2 <- FindClusters(hupsa2, resolution = 0.1, cluster.name = "clusters_0.1")
hupsa2 <- RunUMAP(hupsa2,reduction = "pca",dims = 1:pcSelect, reduction.name = "umap",return.model = T)

DimPlot(hupsa2, group.by = "seurat_clusters", label = T,raster = F)
DimPlot(hupsa2, group.by = "cell_type", label = T,raster = F)

saveRDS(hupsa2, "/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/HuPSA_tumoronly(250430).rds")

#############################################################
data <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/epi_rough(241024).rds")
data <- UpdateSeuratObject(data)

### start
anchor <- FindTransferAnchors(
  reference = hupsa2,
  query = data,
  reference.reduction = "pca",
  normalization.method = "LogNormalize",
  dims = 1:10
)

gc()

data=MapQuery(
  anchorset = anchor,
  query = data,
  reference = hupsa2,
  reference.reduction = "pca",
  refdata = list(
    seurat_clusters = "seurat_clusters",
    cell_type = "cell_type",
    histo = "histo"
  ),
  reduction.model = "umap"
)



SCpubr::do_DimPlot(sample = hupsa,pt.size = 0.1,raster = F,border.size = 20,legend.position = "none",label = F,label.color = "black")+
  coord_cartesian(xlim = c(-15, 15), ylim = c(-15, 15)) 

SCpubr::do_DimPlot(sample = data,reduction = "ref.umap",pt.size = 0.1,raster = F,border.size = 20,legend.position = "none",label = F,label.color = "black")

SCpubr::do_DimPlot(sample = data, reduction = "ref.umap",pt.size = 0.1,raster = F,border.size = 20,legend.position = "none",label = F,
                   group.by = "predicted.cell_type", label.color = "black", split.by = "Sample_ID")+
  coord_cartesian(xlim = c(-15, 15), ylim = c(-15, 15)) 


DimPlot2(data,reduction = "ref.umap", group.by = "predicted.seurat_clusters",split.by = "Sample_ID",theme = theme_umap_arrows())+
  scale_x_continuous(limits = c(-15, 15)) +
  scale_y_continuous(limits = c(-15, 15))

DimPlot2(data,reduction = "ref.umap", group.by = "predicted.cell_type",split.by = "Sample_ID",theme = theme_umap_arrows())+
  scale_x_continuous(limits = c(-15, 15)) +
  scale_y_continuous(limits = c(-15, 15))

DimPlot2(hupsa2, group.by = "cell_type",theme = theme_umap_arrows())+
  scale_x_continuous(limits = c(-15, 15)) +
  scale_y_continuous(limits = c(-15, 15))

DimPlot2(hupsa2, group.by = "cell_type",theme = theme_umap_arrows())

DimPlot(data,reduction = "ref.umap", group.by = "epi_rough")
DimPlot(data,reduction = "ref.umap", split.by = "Sample_ID")

table(hupsa2$histo)



Idents(hupsa) <- "cell_type"
AR_1 <- subset(hupsa, idents = "AdPCa_AR+_1")

DimPlot(AR_1, group.by = "cell_type")+
  coord_cartesian(xlim = c(-15, 15), ylim = c(-15, 15)) 


DimPlot(data,reduction = "ref.umap", group.by = "predicted.cell_type")+
  coord_cartesian(xlim = c(-15, 15), ylim = c(-15, 15)) 

DimPlot(data,reduction = "ref.umap", group.by = "predicted.seurat_clusters",split.by = "Sample_ID")+
  coord_cartesian(xlim = c(-15, 15), ylim = c(-15, 15)) 

DimPlot(hupsa, group.by = "cell_type")+
  coord_cartesian(xlim = c(-15, 15), ylim = c(-15, 15)) 

DimPlot(data,reduction = "ref.umap", group.by = "predicted.seurat_clusters",split.by = "Sample_ID")+
  coord_cartesian(xlim = c(-15, 15), ylim = c(-15, 15)) 









#############################################
meta <- FetchData(hupsa, vars = c("sample", "seurat_clusters"))
colnames(meta) <- c("Sample_ID", "cluster")

meta_percentages <- meta %>%
  group_by(Sample_ID, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%  # Count occurrences
  group_by(Sample_ID) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()


# Determine the order of cell types
# unique_cell_types <- factor(c(0:22))

# Set the levels of cell_type in the meta_percentages data frame
# meta_percentages$cluster <- factor(meta_percentages$cluster, levels = unique_cell_types)

colors <- c("#636363", "#d8b365", "#5ab4ac", "#d73027", "#4575b4",
            "#fb9a99", "#fdbf6f", "#cab2d6", "#b15928", "#87CEFA",
            "#b5ae66", "#ce1256", "#8c510a", "#bf812d", "#35978f",
            "#1a9850", "#313695", "#b2182b", "#fde0ef", "#f7f7f7",
            "#de77ae", "#91bfdb", "#fee08b", "#d53e4f")
ggplot(meta_percentages, aes(x = Sample_ID, y = percentage, fill = cluster)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  labs(y = "Percentage (%)", x = "Sample_ID", title = "Cell Type Proportions by ID") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))+
  coord_flip()













##############
# 3rd
Idents(hupsa2) <- "seurat_clusters"
to_remove <- c("5","7")  
hupsa3 <- subset(hupsa2, idents = to_remove, invert = TRUE)
DimPlot(hupsa3, raster = F)
length(unique(hupsa3$sample))

set.seed(123)
pcSelect <- 10
hupsa3 <- FindNeighbors(hupsa3, reduction = "pca", dims = 1:pcSelect)
hupsa3 <- FindClusters(hupsa3, resolution = 0.1, cluster.name = "clusters_0.1")
hupsa3 <- RunUMAP(hupsa3,reduction = "pca",dims = 1:pcSelect, reduction.name = "umap",return.model = T)

DimPlot(hupsa3, group.by = "seurat_clusters", label = T,raster = F)
DimPlot(hupsa3, group.by = "cell_type", label = F,raster = F)

#-------------------
data <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/epi_rough(241024).rds")
data <- UpdateSeuratObject(data)

### start
anchor <- FindTransferAnchors(
  reference = hupsa3,
  query = data,
  reference.reduction = "pca",
  normalization.method = "LogNormalize",
  dims = 1:10
)

gc()

data=MapQuery(
  anchorset = anchor,
  query = data,
  reference = hupsa3,
  reference.reduction = "pca",
  refdata = list(
    cell_type = "cell_type"
  ),
  reduction.model = "umap"
)

#--------------
colors <- c("#d53e4f", "#bf812d", "#5ab4ac", "#4575b4")

DimPlot(data,reduction = "ref.umap", group.by = "predicted.cell_type")+
  coord_cartesian(xlim = c(-13, 12), ylim = c(-11, 15)) 
DimPlot(hupsa3,reduction = "umap", group.by = "cell_type")+
  coord_cartesian(xlim = c(-13, 12), ylim = c(-11, 15))
DimPlot(hupsa3,reduction = "umap", group.by = "histo")+
  coord_cartesian(xlim = c(-15, 15), ylim = c(-15, 15))

DimPlot(data,reduction = "ref.umap", group.by = "epi_rough", cols = colors)+
  coord_cartesian(xlim = c(-13, 12), ylim = c(-11, 15))

# colors <- c("#91bfdb", "#bf812d", "#d53e4f")
DimPlot(data,reduction = "ref.umap", group.by  = "Sample_ID",pt.size = 0.5,alpha = 0.5)+
  coord_cartesian(xlim = c(-13, 12), ylim = c(-11, 15))

DimPlot(data,reduction = "ref.umap", group.by  = "Sample_ID")


#-----
umap_df <- as.data.frame(Embeddings(data, "ref.umap"))
umap_df$cluster <- data$Sample_ID


values = c("c0"="#d53e4f","c1"="#bf812d","c2"="#5ab4ac","c3"="#4575b4")
values = c("sgNT"="#91bfdb","sgNLN"="#bf812d","sgNLN_WNT"="#d53e4f")

ggplot(umap_df, aes(x = refUMAP_1, y = refUMAP_2, color = cluster)) +
  geom_jitter(width = 0.3, height = 0.3, size = 0.8, alpha = 0.6) +
  scale_color_manual(values = c("sgNT"="#91bfdb","sgNLN"="#bf812d","sgNLN_WNT"="#d53e4f")) +  # change to your actual cluster names
  labs(title = "UMAP with Jittering") +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

