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
library(SeuratExtend)

options(Seurat.object.assay.version = 'v3')

# tumor.epi <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/primary_harmony_p20r0.1(241023).rds")
# DimPlot(tumor.epi, label = T)
tumor.epi <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/epi(250429).rds")
tumor.epi <- RenameIdents(object = tumor.epi,
                          "0" = "c0",
                          "1" = "c1",
                          "2" = "c2",
                          "3" = "c3")

tumor.epi$epi_rough <- Idents(tumor.epi)
Idents(tumor.epi) <- "epi_rough"
DimPlot(tumor.epi, label = T)

saveRDS(tumor.epi, "/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/epi_rough(250429).rds")


#-----------------------------
tumor.epi <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/epi_rough(250429).rds")

colors <- c("c0"="#d53e4f", "c1"="#bf812d", "c2"="#5ab4ac", "c3"="#4575b4")
DimPlot2(tumor.epi, theme = theme_umap_arrows(), label = F, cols = colors)

DimPlot(object = tumor.epi, reduction = "harmony", group.by = "epi_rough", cols = colors)
FeaturePlot(tumor.epi, features = "KEGG_WNT_SIGNALING_PATHWAY1",reduction = "harmony", order = F,cols = c("#FFF1ED", "#3F0A3F"))
FeaturePlot(tumor.epi, features = c("KEGG_WNT_SIGNALING_PATHWAY1"), reduction = "harmony",pt.size = 0.5, label = F, label.color = "white", label.size = 5, order = F)+ 
  scale_color_viridis_c(option = "magma",limits = c(0, 0.25))

gene_list <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/4.1 Add_gene_signature/lineage_240222_aligned.rds")
genesig <- gene_list[c("Jak.genes", "Jak.genes.long", "SOX2.targe.genes", "Lineage.plastic.sig.Nelson")]
genesig<- gene_list$Stem.ACS.signature
# DNMT1, HELLS, LGR6, PTTG1, CDK6, CD44, CD55, KIT, KLF4, NANOG, NES, NOTCH4, OCT4, PDPN, PROM1, SOX2



for(gn in names(gene_list)){
  tumor.epi <- AddModuleScore(tumor.epi, 
                                      features = list(gene_list[[gn]]), 
                                      assay = "RNA",
                                      name = gn)
}
DotPlot(object = tumor.epi, features = paste(names(gene_list), "1", sep=""), scale.by = "size")+
  scale_colour_gradient2(high = "#3F0A3F", mid = "#B6CDE3", low = "white", midpoint = 0)+ 
  RotatedAxis()

#-----------------
colors <- c("#8c510a", "#fdbf6f","#87CEFA","#d8b365","#5ab4ac","#d73027", "#4575b4", 
            "#fb9a99","#636363", "#cab2d6","#b15928","#87CEFA","#b5ae66","#ce1256")  

colors <- c("#636363", "#d8b365", "#5ab4ac", "#d73027", "#4575b4",
            "#fb9a99", "#fdbf6f", "#cab2d6", "#b15928", "#87CEFA",
            "#b5ae66", "#ce1256", "#8c510a", "#bf812d", "#35978f",
            "#1a9850", "#313695", "#b2182b", "#fde0ef", "#f7f7f7",
            "#de77ae", "#91bfdb", "#fee08b", "#d53e4f")

# "#ffff99", "#33a02c"

colors <- c("#8c510a", "#d53e4f","#4575b4","#d8b365","#5ab4ac","#d73027", "#4575b4", 
            "#fb9a99","#636363", "#cab2d6","#b15928","#87CEFA","#b5ae66","#ce1256")  

colors <- c("#d53e4f", "#bf812d", "#5ab4ac", "#4575b4")

DimPlot(tumor.epi, label = F, raster = F, cols = colors)
DimPlot(tumor.epi, reduction = "umap", label = F)

# Set your desired order of samples
desired_order <- c("sgNT", "sgNLN", "sgNLN_WNT")
tumor.epi$Sample_ID <- factor(tumor.epi$Sample_ID, levels = desired_order)
DimPlot(tumor.epi, reduction = "umap", label = FALSE, raster = FALSE, cols = colors, split.by = "Sample_ID")

colors <- c("#91bfdb", "#d8b365", "#B00015")
DimPlot(tumor.epi, reduction = "umap", label = FALSE, raster = FALSE, cols = colors, order = F,pt.size = 0.1,group.by = "Sample_ID")

table(tumor.epi$Patient_ID)


FeaturePlot(tumor.CD4, features = c("CD40LG","IL7R"), label = F, order = T)
FeaturePlot(tumor.CD4, features = c("CCL20","KLRB1","RORA"), label = F, order = T, split.by = "grouping")

FeaturePlot(object = tumor.CD4, features = marker_genes, order = T, label = F, cols = c("#FFF1ED", "#3F0A3F"))

DotPlot(object = tumor.CD4, features = c("CD40LG","IL7R","CD69"), cols = c("#FFF1ED", "#B00015"), scale.by = "size",group.by = "CD4_rough")+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#--------------------
tumor.CD4$CD4_rough <- factor(tumor.CD4$CD4_rough, 
                              levels = c("Tn_c01","Tem_c02","Tcm_c03","Th17_c04","CTL_c05","Tfh_c06","Treg_c07","GIMAP7+T_c08"))

marker_genes <- c("CD3E","CD3D","CD3G","CD4","SELL","LEF1","CCR7","TCF7","FOS","JUN","GADD45B",
                  "CD69","GPR183","CD44","ANXA1","LMNA","CCL20","KLRB1","RORA","GZMK","GZMA","CCL5","NMB","CD200",
                  "TIGIT","IL2RA","FOXP3","TNFRSF9","GIMAP7","TXNIP","TRAF3IP3")
marker_genes <- c("PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","LAYN","CXCL13","TNFRSF9","EOMES","ENTPD1","CD244","CD160","CXCR6","TCF7","TOX","TOX2",
                  "BTLA","NR4A2","IKZF3","BHLHE40","TMIGD2","BATF","ETV1","RUNX3","IRF1","ID2","PRDM1","TRIB1","CSF1","CCL3","CCL4","CCL5",
                  "NFATC1","NFAT5","RGS16")

DotPlot(object = tumor.CD4, features = marker_genes, dot.scale = 6, scale = T,group.by = "grouping")+
  scale_colour_gradient2(high = "#CB2924", mid = "#EBF3E0", low = "#4677AB", midpoint = 0)+
  scale_size(range = c(0.5,6))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
#CB2924 #EBF3E0 #4677AB
# cols = c("#FFF1ED", "#B00015")
# high = "#671215", mid = "#EC674A", low = "#FEF1EB"

#----------------
tumor.CD4_sub <- subset(tumor.CD4, idents = c("Treg_c07"))
marker_genes <- c("IL2RA","IL2RB","IL2RG","EBI3","HAVCR2","LAG3","TIGIT","CTLA4","ENTPD1","LAYN","ICOS","TNFRSF4",
                  "TNFRSF18","TNFRSF9")
DotPlot(object = tumor.CD4, features = marker_genes, dot.scale = 6, scale = T,group.by = "grouping")+
  scale_colour_gradient2(high = "#CB2924", mid = "#EBF3E0", low = "#4677AB", midpoint = 0)+
  scale_size(range = c(0.5,6))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#-----------------
# find marker genes
DefaultAssay(tumor.epi) <- "RNA"
sc.markers <- FindAllMarkers(object = tumor.epi,
                             only.pos = F,
                             min.pct = 0.1,
                             logfc.threshold = 0.01)

# Reorder columns and sort by log2 fold change   
ann_markers <- sc.markers[ , c(6, 7, 2:5, 1)]

ann_markers <- ann_markers %>%
  dplyr::arrange(dplyr::desc(abs(avg_log2FC))) 

write.table(ann_markers, "/work/SCCC/s228508/5.2.1 CD4/CD4_rough(240818)_markergenes.csv", sep = ",")

#------------------
gene_list <- readRDS("/work/SCCC/s228508/4.1 Add_gene_signature/CD4_sig_aligned(240609).rds")
gene_list <- readRDS("/work/SCCC/s228508/choushi/Endocytosis.rds")

kegg <- readRDS("/work/SCCC/s228508/7.0 Myeloid/KEGG(240929).rds")
for (i in names(kegg)) {
  kegg[[i]] <- unique(kegg[[i]])
}
gene_list <- kegg$KEGG_WNT_SIGNALING_PATHWAY

tumor.epi <- AddModuleScore(tumor.epi, 
                            features = list(gene_list), 
                            assay = "RNA",
                            name = "KEGG_WNT_SIGNALING_PATHWAY")
DotPlot(object = tumor.epi, features = "KEGG_WNT_SIGNALING_PATHWAY1", scale = T)+
  scale_colour_gradient2(high = "#3F0A3F", mid = "#B6CDE3", low = "white", midpoint = 0)+ 
  RotatedAxis()
VlnPlot(object = tumor.epi, features = "KEGG_WNT_SIGNALING_PATHWAY1")
FeaturePlot(object = tumor.epi, features = "KEGG_WNT_SIGNALING_PATHWAY1")


for(gn in names(gene_list)){
  tumor.epi <- AddModuleScore(tumor.epi, 
                              features = list(gene_list[[gn]]), 
                              assay = "RNA",
                              name = gn)
}

DotPlot(object = tumor.epi, features = paste(names(gene_list), "1", sep=""), scale = T)+
  scale_colour_gradient2(high = "#3F0A3F", mid = "#B6CDE3", low = "white", midpoint = 0)+ 
  RotatedAxis()
#------------------
gene_list <- readRDS("/work/SCCC/s228508/4.1 Add_gene_signature/CD4_function_sig_aligned(240818).rds")
for(gn in names(gene_list)){
  tumor.CD4 <- AddModuleScore(tumor.CD4, 
                              features = list(gene_list[[gn]]), 
                              assay = "RNA",
                              name = gn)
}

DotPlot(object = tumor.CD4, features = paste(names(gene_list), "1", sep=""), scale = T)+
  scale_colour_gradient2(high = "#671215", mid = "#EC674A", low = "#FEF1EB", midpoint = 1)+ 
  RotatedAxis()


#--------------
immune_sig <- readRDS("/work/SCCC/s228508/4.1 Add_gene_signature/immune_sig_aligned(20240330).rds")

selected_indices <- c("Terminal exhaustion","Exhaustion","TIL dysfunction","Effector","Progenitor exhausted","CD4 T cell dysfunction in melanoma")
immune_sig2 <- immune_sig[selected_indices]

tumor.CD4_sub <- subset(tumor.CD4, idents = "CTL_c05")
for(gn in names(immune_sig2)){
  tumor.CD4_sub <- AddModuleScore(tumor.CD4_sub, 
                                  features = list(immune_sig2[[gn]]), 
                                  assay = "RNA",
                                  name = gn)
}
DotPlot(object = tumor.CD4_sub, features = paste(names(immune_sig2), "1", sep=""), 
        cols = c("#FEF1EB", "#3F0A3F"), scale.by = "size",group.by = "grouping")+RotatedAxis()

#------------

# tumor.CD8 <- subset(tumor.CD8, idents = c("GZMK+ T","Terminal Exh"))
for(gn in names(effector_lsit)){
  tumor.CD4 <- AddModuleScore(tumor.CD4, 
                              features = list(effector_lsit[[gn]]), 
                              assay = "RNA",
                              name = gn)
}

DotPlot(object = tumor.CD4, features = paste(names(effector_lsit), "1", sep=""), dot.scale = 6, scale = T,group.by = "grouping")+
  scale_colour_gradient2(high = "#3F0A3F", mid = "#B6CDE3", low = "white", midpoint = 0)+ 
  scale_size(range = c(0.5,6))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))




############################################
# heatmap
#------------------
canonical_WNT <- c("HDAC11","PPARD","KREMEN2","PRKCZ","CUL3","RHOA","KLHL12","FLNA","NOTCH1","NCSTN","HDAC5","FZD4","GNAI1","PSEN2","LRP5",
                   "ATP6AP2","WNT5A","WNT3A","CTNNB1","PORCN","FZD1","FZD5","JAG1","JAG2","LEF1","TCF7","LGR4","TCF7L2","ATOH1","KLF4")
canonical_AR <- c("SLC38A2","PMEPA1","ZBTB10","KLK2","KLK3","GUCY1A1","C1orf116","HERC3","IQGAP2","TMPRSS2","ABHD2","NKX3-1","ELL2","ABCC4",
                  "ELOVL5","HOMER2")

gene_cell_exp <- AverageExpression(tumor.epi,
                                   features = canonical_AR,
                                   group.by = "epi_rough",
                                   slot = "data") 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
marker_exp <- scale(t(gene_cell_exp),scale = T,center = T) # scale() scales data according to the columns!! 
marker_exp <- as.matrix(marker_exp)
# marker_exp <- marker_exp[c("pDC","cDC1","cDC2","LAMP3+ DC","Monocytes","TAM_c07_THBS1","TAM_c01_TM4SF1","TAM_c02_APOD","TAM_c03_transitory",
#                            "TAM_c04_FOLR2","TAM_c05_IL1B","TAM_c06_TREM2","TAM_c08_SPP1"),]

color_mapping <- colorRamp2(c(-2, 0, 2), c("#4575B4","white","#D62F25"))

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



genes.zscore <- CalcStats(
  tumor.epi,
  features = canonical_AR,
  group.by = "epi_rough",
  order = "p")

Heatmap(genes.zscore, lab_fill = "zscore")


#####################################################
## cell proportion
#------------
# tumor.CD4[["Idents"]] <- Idents(tumor.CD4)

meta <- FetchData(tumor.epi, vars = c("Sample_ID", "epi_rough"))
colnames(meta) <- c("Sample_ID", "cluster")
meta$Sample_ID <- factor(meta$Sample_ID, levels = c("sgNLN_WNT","sgNLN","sgNT"))

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

# colors <- c("#b15928", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#5ab4ac")

ggplot(meta_percentages, aes(x = Sample_ID, y = percentage, fill = cluster)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  labs(y = "Percentage (%)", x = "Sample_ID", title = "Cell Type Proportions by ID") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10)) +
  coord_flip()

# three big bars
#------------
meta_percentages <- meta %>% 
  count(Patient_ID, Group_Name,cluster) %>%  
  group_by(Group_Name) %>%          
  mutate(percentage = n/sum(n)) %>%  
  ungroup()
# meta_percentages$percentage <- log2(meta_percentages$percentage)

# colors <- c("#1f78b4","#5ab4ac","#6a3d9a","#B07AA1","#b2df8a","#60CEF4","#E15759","#b15928","#fb9a99","#fdbf6f")
ggplot(meta_percentages, aes(x = Group_Name, y = percentage, fill = cluster)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  labs(y = "Percentage", x = "grouping", title = "Cell type proportions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#----------
# separate bars
table(tumor.CD4$Patient_ID)
meta2 <- meta[!(meta$Patient_ID %in% c("ADT_10879_F","ADT_31230_F","ADT_31666_F",
                                       "ADT_31953_F","ADT_32032_F","ADT_32126_F","ADT_32360_Fcol",
                                       "ADT_32378_F","ADT_32413_F","ADT_32667_F","ADT_33005_Fcol")), ]

meta_percentages <- meta2 %>% 
  count(Patient_ID, Group_Name,cluster) %>%  
  group_by(Patient_ID) %>%          
  mutate(percentage = n/sum(n)) %>%  
  ungroup()
# meta_percentages$percentage <- log2(meta_percentages$percentage)

# "Tn_c01","Tem_c02","Tcm_c03","Th17_c04","CTL_c05","Tfh_c06","Treg_c07","GIMAP7+T_c08"
colors <- c('#f6bd60', '#f28482', "#84a59d")#, "#B07AA1")
select <- subset(meta_percentages, meta_percentages$cluster == "GIMAP7+T_c08")
ggplot(select, aes(x = Group_Name, y = percentage, fill = Group_Name)) +
  geom_boxplot() +
  geom_dotplot(aes(fill = Group_Name), binaxis = "y", stackdir = "center", dotsize = 0.8) +
  labs(y = "Fraction in CD4 T", title = "GIMAP7+T_c08") +
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#--------------------------------
# find marker genes
DefaultAssay(tumor.CD4) <- "RNA"
Idents(tumor.CD4) <- "grouping"

ribosomal_genes <- grep("^RPS|^RPL", rownames(tumor.CD4), value = TRUE)
seurat_obj <- subset(tumor.CD4, features = setdiff(rownames(tumor.CD4), ribosomal_genes))

sc.markers <- FindAllMarkers(object = seurat_obj,
                             only.pos = F,
                             min.pct = 0.2,
                             logfc.threshold = 0.01)

ann_markers <- sc.markers[ , c(6, 7, 2:5, 1)]

ann_markers <- ann_markers %>%
  dplyr::arrange(dplyr::desc(abs(avg_log2FC))) 

#-------
ribosomal_genes <- grep("^RPS|^RPL", rownames(tumor.CD4), value = TRUE)
seurat_obj <- subset(tumor.CD4, features = setdiff(rownames(tumor.CD4), ribosomal_genes))

sc.markers <- FindMarkers(object = seurat_obj,
                          ident.1 = "PCa_AR_high",
                          ident.2 = "ADT_lineage_high",
                          only.pos = F,
                          min.pct = 0.1,
                          logfc.threshold = 0.01)
# Reorder columns and sort by log2 fold change   
ann_markers2 <- sc.markers[ , c(6, 7, 2:5, 1)]

ann_markers2 <- ann_markers2 %>%
  dplyr::arrange(dplyr::desc(abs(avg_log2FC))) 

write.table(ann_markers, "/work/SCCC/s228508/5.2 TNK/marker_CD8_harmony_G(240405)_2.csv", sep = ",")








#####################################################
# Cell cycle
cc.genes <- read.table("/gpfs/gibbs/project/mu/qx75/Grey/1.0 integration/cc_genes.csv", header=T, sep=",", check.names=F)
remove_empty <- function(x) x[x != ""]
cc.genes <- lapply(cc.genes, remove_empty)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

tumor.epi <- CellCycleScoring(tumor.epi, g2m.features = g2m.genes, s.features = s.genes)

# tumor.epi <- RunTSNE(tumor.epi, reduction = "harmony", dims = 1:pcSelect)

desired_order <- c("sgNT", "sgNLN", "sgNLN_WNT")
tumor.epi$Sample_ID <- factor(tumor.epi$Sample_ID, levels = desired_order)

colors <- c("#91bfdb", "#d8b365", "#636363")
DimPlot2(tumor.epi, theme = theme_umap_arrows(), group.by = "Phase",label = F, cols = colors)


#----------
# cell cycle proportions
meta <- FetchData(tumor.epi, vars = c("Sample_ID", "Phase"))
colnames(meta) <- c("Sample_ID", "cluster")

colors <- c("#91bfdb", "#d8b365", "#636363")

meta_percentages <- meta %>%
  group_by(Sample_ID, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%  # Count occurrences
  group_by(Sample_ID) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()

# meta_percentages$percentage <- log2(meta_percentages$percentage)

# colors <- c("#1f78b4","#5ab4ac","#6a3d9a","#B07AA1","#b2df8a","#60CEF4","#E15759","#b15928","#fb9a99","#fdbf6f")
ggplot(meta_percentages, aes(x = Sample_ID, y = percentage, fill = cluster)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  labs(y = "Percentage", x = "grouping", title = "Cell type proportions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_flip()











##########################################################################################
kegg_df <- msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
kegg <- kegg_df %>% split(x = .$gene_symbol, f = .$gs_name)

tumor.epi <- AddModuleScore(tumor.epi, 
                             features = list(kegg[["KEGG_WNT_SIGNALING_PATHWAY"]]), 
                             assay = "RNA",
                             name = "KEGG_WNT_SIGNALING_PATHWAY")

#-------------------------------
meta <- FetchData(tumor.epi, vars = c("KEGG_WNT_SIGNALING_PATHWAY1","Sample_ID"))
colnames(meta) <- c("KEGG_WNT_SIGNALING_PATHWAY1","Sample_ID")

# write.table(meta, "/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/wnt_table2.csv", sep = ",")
# wilcox.test(REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES1 ~ SYNCRIP_patient_class, data = meta)

# VlnPlot(object = seurat_epi, features = "REACTOME_CYTOSOLIC_SENSORS_OF_PATHOGEN_ASSOCIATED_DNA1", pt.size= 0, group.by = "SYNCRIP_patient_class",
#         cols = c("SYNCRIP_high"="#D6606B", "SYNCRIP_low"="#4575b4"))+
#   geom_violin(trim = FALSE, scale = "width") + 
#   # geom_jitter(position = position_jitter(width = 0.2), size = 0.001, alpha = 0.8) +
#   theme_classic() +
#   theme(
#     axis.title.x = element_text(size = 12),
#     axis.title.y = element_text(size = 12),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     plot.title = element_text(size = 14, face = "bold")
#   ) +
#   xlab("Genotype") +
#   ylab("Gene Expression")
colors <- c("#d53e4f", "#bf812d", "#5ab4ac", "#4575b4")
colors <- c("#91bfdb", "#bf812d", "#d53e4f")

VlnPlot(object = tumor.epi, 
        features = "KEGG_WNT_SIGNALING_PATHWAY1", 
        pt.size = 0, group.by = "Sample_ID",
        cols = colors)+
  geom_violin(trim = FALSE, scale = "width") +  # Violin plot
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white", color = "black") +  # Boxplot inside
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


#------------
HALLMARK_WNT_BETA_CATENIN_SIGNALING
GOBP_CANONICAL_WNT_SIGNALING_PATHWAY
KEGG_WNT_SIGNALING_PATHWAY
PID_WNT_SIGNALING_PATHWAY
HALLMARK_IL6_JAK_STAT3_SIGNALING

kegg_df <- msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
kegg_df <- msigdbr(species = "human", category = "H") 
kegg_df <- msigdbr(species = "human", category = "C5", subcategory = "GO:BP")
kegg_df <- msigdbr(species = "human", category = "C2", subcategory = "CP:PID")

kegg <- kegg_df %>% split(x = .$gene_symbol, f = .$gs_name)

tumor.epi <- AddModuleScore(tumor.epi, 
                            features = list(unique(kegg[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]])), 
                            assay = "RNA",
                            name = "HALLMARK_IL6_JAK_STAT3_SIGNALING")



# Specifying genes and cells of interest
genes <- c("CD3D", "CD14", "CD79A")
cells <- WhichCells(pbmc, idents = c("B cell", "CD8 T cell", "Mono CD14"))

desired_order <- c("sgNT", "sgNLN", "sgNLN_WNT")
tumor.epi$Sample_ID <- factor(tumor.epi$Sample_ID, levels = desired_order)
colors <- c("#d53e4f", "#bf812d", "#5ab4ac", "#4575b4")
colors <- c("#91bfdb", "#d8b365", "#636363")
# Violin plot with statistical analysis
VlnPlot2(
  tumor.epi,
  features = "HALLMARK_IL6_JAK_STAT3_SIGNALING1",
  group.by = "epi_rough",
  # cells = cells,
  cols = colors,
  stat.method = "wilcox.test")
