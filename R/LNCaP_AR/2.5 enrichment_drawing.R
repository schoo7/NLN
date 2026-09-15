rm(list=ls())
library(RColorBrewer)
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(colorspace)  
library(ggplot2)
library(dplyr)
library(harmony)
library(cowplot)
library(ComplexHeatmap)
library(circlize)

# tumor.epi <- readRDS("/work/SCCC/s228508/choushi/tumor.epi_top1000_p20_res0.1.rds")
# colors <- c("0"='#84a59d', "1"='#f28482', "2"="#f6bd60", "3"="#B07AA1")
# 
# pdf("/work/SCCC/s228508/choushi/figs/UMAP_tumor_only.pdf", width=5, height=4)
# DimPlot(tumor.epi, label = F, raster = F, cols = colors)
# dev.off()
# 
# pdf("/work/SCCC/s228508/choushi/figs/dotplot_NLN.pdf", width=3.5, height=3.5)
# DotPlot(tumor.epi, features = "NLN")+
#   scale_colour_gradient2(high = "#3F0A3F", mid = "#B6CDE3", low = "white", midpoint = 0)+ 
#   RotatedAxis()
# dev.off()

# tumor.epi <- RenameIdents(object = tumor.epi, 
#                           "0" = "NLN_low",
#                           "1" = "NLN_low",
#                           "2" = "NLN_high",
#                           "3" = "NLN_low")
# 
# tumor.epi$NLN_class <- Idents(tumor.epi)
# Idents(tumor.epi) <- "NLN_class"
# saveRDS(tumor.epi, "/work/SCCC/s228508/choushi/tumor.epi_NLN_class(240511).rds")

tumor.epi <- readRDS("/work/SCCC/s228508/grey_Mu_lab/choushi/primary_harmony(240929).rds")

colors <- c("NLN_low"="#D6606B", "NLN_high"="#4575b4")
DimPlot(tumor.epi, label = T, raster = F,cols = colors, split.by = "Patient_ID") 
# pdf("/work/SCCC/s228508/choushi/figs/UMAP.pdf", width=5, height=4)
DimPlot(tumor.epi, label = T, raster = F) 
# dev.off()

DefaultAssay(tumor.epi) <- "RNA"
showGenes = c("NLN") 
FeaturePlot(object = tumor.epi, features = showGenes, pt.size = 1,order = T, label = , cols = c("#FFF1ED", "#3F0A3F"))
VlnPlot(object = tumor.epi, features = "NLN",pt.size =0)
DotPlot(object = tumor.epi, features = showGenes, cols = c("#FFF1ED", "#B00015"), col.max = 2.5, scale.by = "size")+ 
  scale_size(range = c(0.5,6))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#----------------
#load HAKKMARK, kegg
m_df<- msigdbr(species = "human",  category = "H" )
geneSets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

kegg <- readRDS("/work/SCCC/s228508/7.0 Myeloid/KEGG(240929).rds")
for (i in names(kegg)) {
  kegg[[i]] <- unique(kegg[[i]])
}

for(gn in names(kegg)){
  tumor.epi <- AddModuleScore(tumor.epi, 
                              features = list(kegg$KEGG_WNT_SIGNALING_PATHWAY), 
                              assay = "RNA",
                              name = "KEGG_WNT_SIGNALING_PATHWAY")
}

meta <- FetchData(tumor.epi, vars = c("HALLMARK_APOPTOSIS1","NLN_class"))
colnames(meta) <- c("HALLMARK_APOPTOSIS1","NLN_class")
wilcox.test(HALLMARK_APOPTOSIS1 ~ NLN_class, data = meta, paired=FALSE)


VlnPlot(object = tumor.epi, features = "KEGG_WNT_SIGNALING_PATHWAY1", pt.size= 0)+
  geom_violin(trim = FALSE, scale = "width") + 
  geom_jitter(position = position_jitter(width = 0.2), size = 0.001, alpha = 0.8) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  xlab("Genotype") +
  ylab("Gene Expression")
ggsave("/work/SCCC/s228508/choushi/figs/apoptosis_1.pdf", width = 7, height = 6)

VlnPlot(object = tumor.epi, features = "HALLMARK_APOPTOSIS1", pt.size= 0)+
  geom_boxplot(width= 0.2,outlier.size=0,fill="white")
ggsave("/work/SCCC/s228508/choushi/figs/apoptosis_2.pdf", width = 6, height = 6)

#-----------------
gene_list <- readRDS("/work/SCCC/s228508/choushi/gene_list_wnt.rds")

for(gn in names(gene_list)){
  tumor.epi <- AddModuleScore(tumor.epi, 
                              features = list(gene_list[[gn]]), 
                              assay = "RNA",
                              name = gn)
}

meta <- FetchData(tumor.epi, vars = c("GOBP_CLATHRIN_DEPENDENT_ENDOCYTOSIS1", "HALLMARK_WNT_BETA_CATENIN_SIGNALING1","NLN_class"))
colnames(meta) <- c("GOBP_CLATHRIN_DEPENDENT_ENDOCYTOSIS1", "HALLMARK_WNT_BETA_CATENIN_SIGNALING1","NLN_class")
wilcox.test(HALLMARK_WNT_BETA_CATENIN_SIGNALING1 ~ NLN_class, data = meta, paired=FALSE)

DotPlot(object = tumor.epi, features = paste(names(gene_list), "1", sep=""), cols = c("#FFF1ED", "#B00015"), scale.by = "size")+
  scale_size(range = c(0.5,6))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#-------
for(gn in names(gene_list)){
  VlnPlot(object = tumor.epi, features = paste(gn, "1", sep=""), pt.size= 0, cols = c("NLN_low"="#D6606B", "NLN_high"="#4575b4"))+
    geom_violin(trim = FALSE, scale = "width") + 
    geom_jitter(position = position_jitter(width = 0.2), size = 0.001, alpha = 0.8) +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    xlab("Genotype") +
    ylab("Gene Expression")
  ggsave(paste0("/work/SCCC/s228508/choushi/figs/", gn,".pdf"), width = 7, height = 6)
}

#--------
for(gn in names(gene_list)){
  VlnPlot(object = tumor.epi, features = paste(gn, "1", sep=""), pt.size= 0, cols = c("NLN_low"="#D6606B", "NLN_high"="#4575b4"))+
    geom_boxplot(width= 0.2,outlier.size=0,fill="white")
  ggsave(paste0("/work/SCCC/s228508/choushi/figs/", gn,".pdf"), width = 6, height = 6)
}


for(gn in names(gene_list)){
  FeaturePlot(object = tumor.epi, features = paste(gn, "1", sep=""), label = TRUE, repel = TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  ggsave(paste0("/work/SCCC/s228508/choushi/figs/", gn,".pdf"),width=5, height = 5)
}


for(gn in names(gene_list)){
  range_color = range(tumor.epi[[paste(gn, "1", sep="")]])
  FeaturePlot(object=tumor.epi, features=paste(gn, "1", sep=""), 
              label = TRUE, raster = T, keep.scale = "all", order=T) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n=11, name="RdBu")),
                           limits=range_color)&
    theme(legend.position = "right")
  ggsave(paste0("/work/SCCC/s228508/choushi/figs/",gn,".pdf"),width=5, height = 5)
}

#---------
meta <- FetchData(tumor.epi, vars = c("Patient_ID", "NLN_class","HALLMARK_WNT_BETA_CATENIN_SIGNALING1","GOBP_EXOCYTOSIS1",
                                      "GOBP_CLATHRIN_DEPENDENT_ENDOCYTOSIS1","HALLMARK_MITOTIC_SPINDLE1"))
colnames(meta) <- c("Patient_ID", "NLN_class","HALLMARK_WNT_BETA_CATENIN_SIGNALING1","GOBP_EXOCYTOSIS1",
                    "GOBP_CLATHRIN_DEPENDENT_ENDOCYTOSIS1","HALLMARK_MITOTIC_SPINDLE1")

colors <- c('#f28482', "#84a59d")#, "#B07AA1")
ggplot(meta, aes(x = NLN_class, y = HALLMARK_MITOTIC_SPINDLE1, fill = NLN_class)) +
  geom_boxplot() +
  # geom_dotplot(aes(fill = NLN_class), binaxis = "y", stackdir = "center", dotsize = 0.8) +
  labs(y = "signature score", title = "HALLMARK_MITOTIC_SPINDLE") +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# lineage
#------------------------------------------------
lineage_list <- readRDS("/work/SCCC/s228508/4.1 Add_gene_signature/lineage_240222_aligned.rds")

# lineage_list[["Human.ARPC"]] = NULL
# lineage_list[["Human.NEPC"]] = NULL

# lineage_features <- c()
# for (i in (1:length(lineage_list))) {
#   lineage_features <- append(lineage_features, lineage_list[[i]])
# }
# 
# lineage_features <- unique(lineage_features)

for(gn in names(lineage_list)){
  tumor.epi <- AddModuleScore(tumor.epi, 
                              features = list(lineage_list[[gn]]), 
                              assay = "RNA",
                              name = gn)
}

DotPlot(object = tumor.epi, features = paste(names(lineage_list), "1", sep=""))+ RotatedAxis()




# heatmap
#--------------
wnt_features <- c()
for (i in c(1,3)) {
  wnt_features <- append(wnt_features, gene_list[[i]])
}

wnt_features <- unique(wnt_features)

DotPlot(object = tumor.epi, features = wnt_features, cols = c("#FFF1ED", "#B00015"), scale.by = "size")+
  scale_size(range = c(0.5,6))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#-------------
DefaultAssay(tumor.epi) <- "RNA"
gene_cell_exp <- AverageExpression(tumor.epi,
                                   features = wnt_features,
                                   group.by = "NLN_class",
                                   slot = "data") 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
# marker_exp <- t(as.matrix(gene_cell_exp))
marker_exp <- scale(t(gene_cell_exp),scale = T,center = T) # scale() scales data according to the columns!! 

# numeric_rownames <- as.numeric(rownames(marker_exp))
# marker_exp <- marker_exp[order(numeric_rownames, decreasing = TRUE), ]
# marker_exp <- marker_exp[c(3,2,1),]
color_mapping <- colorRamp2(c(-1, 0, 1), c("#4575B4","white","#D62F25")) #FAF5BD
# color_mapping <- colorRamp2(c(0, 0.5), c("white","#D62F25")) #FAF5BD
#--------
Heatmap(marker_exp,
        cluster_rows = F,
        cluster_columns = F,
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





