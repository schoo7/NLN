rm(list=ls())

library(RColorBrewer)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(pheatmap)
library(homologene)

# cc.genes <- read.table("/work/SCCC/s228508/1.0 integration/cc_genes.csv", header=T, sep=",", check.names=F)
# remove_empty <- function(x) x[x != ""]
# cc.genes <- lapply(cc.genes, remove_empty)
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes

options(Seurat.object.assay.version = 'v3')
# seurat_list <- readRDS("/work/SCCC/s228508/grey_Mu_lab/choushi/diet_seurat_list_all_genes(240928).rds")

seurat_list <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/diet_seurat_list(241021).rds")
seurat.primary <- merge(x = seurat_list[[1]], y = seurat_list[-1])

# seurat_list <- lapply(X = seurat_list, FUN = function(x){
#   x <- NormalizeData(x)
#   # x <- CellCycleScoring(x, g2m.features = g2m.genes, s.features = s.genes)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# }) 
# 
# features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 2000, verbose = TRUE)

# genes_to_remove <- grep("^LINC", rownames(seurat.primary), value = TRUE)
# seurat.primary <- subset(seurat.primary, features = setdiff(rownames(seurat.primary), genes_to_remove))

#-------------------
DefaultAssay(seurat.primary) <- "RNA"
seurat.primary <- NormalizeData(seurat.primary)
# seurat.primary <- CellCycleScoring(seurat.primary, g2m.features = g2m.genes, s.features = s.genes)

seurat.primary <- FindVariableFeatures(seurat.primary, nfeatures = 2000)

all.genes <- rownames(seurat.primary)
seurat.primary <- ScaleData(seurat.primary, features = all.genes, verbose = T)#vars.to.regress = c("S.Score","G2M.Score")
seurat.primary <- RunPCA(seurat.primary,npcs = 50, verbose = T)
DimPlot(object = seurat.primary, reduction = "pca", group.by = "Sample_ID")

seurat.harmony <- RunHarmony(seurat.primary, group.by.vars = c("Sample_ID"), plot_convergence = TRUE)
colors <- c("#91bfdb", "#d8b365", "#b2182b")
DimPlot(object = seurat.harmony, reduction = "harmony", group.by = "Sample_ID", cols = colors)
# VlnPlot(object = seurat.harmony, features = "harmony_1", group.by = "Sample_ID")

ElbowPlot(seurat.primary, ndims = 50)

#-----------
set.seed(123)
pcSelect = 20
seurat.harmony <- FindNeighbors(seurat.harmony, reduction = "harmony", dims = 1:pcSelect)
seurat.harmony <- FindClusters(seurat.harmony, resolution = 0.1) 
# seurat.harmony <- RunUMAP(seurat.harmony, reduction = "harmony", dims = 1:pcSelect)

seurat.harmony <- RunTSNE(seurat.harmony, reduction = "harmony", dims = 1:pcSelect)

colors <- c("#d53e4f", "#bf812d", "#5ab4ac", "#4575b4")
TSNEPlot(seurat.harmony, label = F, cols = colors)
DimPlot2(seurat.harmony, theme = theme_umap_arrows(), label = F, cols = colors)
# saveRDS(seurat.harmony, "/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/epi(250429).rds")


desired_order <- c("sgNT", "sgNLN", "sgNLN_WNT")
seurat.harmony$Sample_ID <- factor(seurat.harmony$Sample_ID, levels = desired_order)

colors <- c("#91bfdb", "#d8b365", "#b2182b")
TSNEPlot(seurat.harmony, label = F, group.by = "Sample_ID",cols = colors)

# saveRDS(seurat.primary, "/work/SCCC/s228508/grey_Mu_lab/choushi/primary_harmony_no_integration(241023).rds")

DimPlot(seurat.harmony, reduction = "umap", label = T, raster= F) 

DimPlot(seurat.primary, reduction = "umap", label = T, raster = F, split.by = "Sample_ID")

FeaturePlot(object = seurat.primary, features = "KRT5", split.by = "Sample_ID", label = T, order = T)
DotPlot(object = seurat.primary, features = c("KIF11","NLN","WNT3A"), group.by = "Sample_ID")+ RotatedAxis()

DefaultAssay(seurat.harmony) <- "RNA" 
genes <- c("Atp6v0d2", "Foxi1","Ascl3","Drapc1")
genes <- c("EPCAM","KRT8","KRT18","ASCL1","NEUROD1","CHGA","CHGB")

genes <- c("Col1a1","Col3a1","Col1a2","Col6a1","Col6a2","Cldn5","Flt1","Myh11")
genes <- c("Sox10","Mpz","Mbp","Gfap","Egr2","Mag","S100b","Plp1")

gene_list <- readRDS("/work/SCCC/s228508/4.1 Add_gene_signature/lineage_240222_aligned.rds")
AR <- gene_list$AR.Score.Gene

NE <- c("CHGA","CHGB","ASCL1","ENO2","SYP","SCG3","SCN3A","PCSK1","ELAVL4","NKX2-1","POU3F2", "INSM1", "UCHL1", "NCAM1", "NRXN1", "CNTN4", "KCNB2")
canonical_AR <- c("SLC38A2","SCD","GUCY1A1","ELL2","KLK2","HOMER2","KLK3","NKX3-1","ABHD2","ELOVL5","ACSL3","TMPRSS2","FKBP5","ABCC4","AKAP12","STK39",
                  "ZBTB10","FADS1","CENPN","IQGAP2")
canonical_WNT <- c("HDAC11","PPARD","KREMEN2","PRKCZ","CUL3","RHOA","KLHL12","FLNA","NOTCH1","NCSTN","HDAC5","FZD4","GNAI1","PSEN2","LRP5",
                   "ATP6AP2","WNT5A","WNT3A","CTNNB1","PORCN","FZD1","FZD5","JAG1","JAG2","LEF1","TCF7","LGR4","TCF7L2","ATOH1","KLF4")

EMT <- c("CDH2","VIM","SNAI1","SNAI2","TWIST1","SMAD2","TGFB1","EPAS1")

FeaturePlot(object = seurat.primary, features = EMT, label = T,order = T)

DotPlot(object = seurat.primary, features = canonical_WNT)+
  scale_colour_gradient2(high = "#3F0A3F", mid = "#B6CDE3", low = "white", midpoint = 0)+ 
  RotatedAxis()

#-----------
# heatmap
gene_cell_exp <- AverageExpression(seurat.primary,
                                   features = canonical_WNT,
                                   group.by = "Sample_ID",
                                   slot = "data") 

gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
marker_exp <- scale(t(gene_cell_exp),scale = T,center = T) # scale() scales data according to the columns!! 
# rownames(marker_exp) <- c("c0","c1","c2","c3","c4")

marker_exp <- marker_exp[c("sgNT","sgNLN-WNT","sgNLN","sgNLN-KIF"),]

color_mapping <- colorRamp2(c(-2, 0, 2), c("#4575B4","white","#D62F25"))

#--------
Heatmap(marker_exp,
        cluster_rows = T,
        cluster_columns = T,
        show_column_names = T,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '),
        col = color_mapping,
        border = 'white',
        rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        width = ncol(marker_exp)*unit(5, "mm"), 
        height = nrow(marker_exp)*unit(5, "mm"),
        column_names_rot = 90,
        row_names_rot = 0,
        row_names_side = "left")





#------------------
gene_list <- readRDS("/work/SCCC/s228508/4.1 Add_gene_signature/lineage_240222_aligned.rds")
gene_list <- readRDS("/work/SCCC/s228508/choushi/Endocytosis.rds")
gene_list <- readRDS("/work/SCCC/s228508/7.0 Myeloid/HALLMARK(240816).rds")
gene_list <- gene_list$HALLMARK_WNT_BETA_CATENIN_SIGNALING

kegg <- readRDS("/work/SCCC/s228508/7.0 Myeloid/KEGG(240929).rds")
for (i in names(kegg)) {
  kegg[[i]] <- unique(kegg[[i]])
}
gene_list <- kegg$KEGG_WNT_SIGNALING_PATHWAY

for(gn in names(gene_list)){
  seurat.primary <- AddModuleScore(seurat.primary, 
                                   features = list(gene_list[[gn]]), 
                                   assay = "RNA",
                                   name = gn)
}
DotPlot(object = seurat.primary, features = paste(names(gene_list), "1", sep=""), scale = T)+
  scale_colour_gradient2(high = "#3F0A3F", mid = "#B6CDE3", low = "white", midpoint = 0)+ 
  RotatedAxis()


seurat.primary <- AddModuleScore(seurat.primary, 
                                 features = list(gene_list), 
                                 assay = "RNA",
                                 name = "KEGG_WNT_SIGNALING_PATHWAY")
DotPlot(object = seurat.primary, features = "KEGG_WNT_SIGNALING_PATHWAY1", scale = T)+
  scale_colour_gradient2(high = "#3F0A3F", mid = "#B6CDE3", low = "white", midpoint = 0)+ 
  RotatedAxis()


#---------------
DefaultAssay(seurat.primary) <- "RNA" 
Idents(seurat.primary) <- "seurat_clusters"

markers <- FindMarkers(seurat.primary, ident.1 = "1", #ident.2 = "sgNT", 
                       group.by = 'seurat_clusters',
                       only.pos = F,
                       min.pct = 0.1,
                       logfc.threshold = 0.01)

# write.table(ann_markers, "/work/SCCC/s228508/grey_Mu_lab/choushi/markergenes.csv", sep = ",")

sc.markers <- FindAllMarkers(object = seurat.primary,
                             only.pos = F,
                             min.pct = 0.1,
                             logfc.threshold = 0.01)

# Reorder columns and sort by log2 fold change   
ann_markers <- sc.markers[ , c(6, 7, 2:5, 1)]

ann_markers <- ann_markers %>%
  dplyr::arrange(dplyr::desc(abs(avg_log2FC))) 

write.table(ann_markers, "/work/SCCC/s228508/grey_Mu_lab/choushi/markergenes.csv", sep = ",")


#------------------------
top3 <- sc.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
DoHeatmap(seurat.primary, top3$gene, size = 3)


