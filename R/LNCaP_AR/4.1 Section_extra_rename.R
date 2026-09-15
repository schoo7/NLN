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

hupsa <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/HuPSA_tumoronly(250430).rds")
hupsa <- RenameIdents(object = hupsa,
                          "0" = "AdPCa-1",
                          "1" = "AdPCa-1",
                          "2" = "NEPC-1",
                          "3" = "AdPCa-3",
                          "4" = "Progenitor-1",
                          "5" = "AdPCa-6",
                          "6" = "AdPCa-2",
                          "7" = "AdPCa-7",
                          "8" = "AdPCa-5",
                          "9" = "NEPC-2",
                          "10" = "AdPCa-4",
                          "11" = "Progenitor-2")

hupsa$epi <- Idents(hupsa)
Idents(hupsa) <- "epi"

desired_order <- c("AdPCa-1","AdPCa-2","AdPCa-3","AdPCa-4","AdPCa-5","AdPCa-6","AdPCa-7","NEPC-1","NEPC-2","Progenitor-1","Progenitor-2")
hupsa$epi <- factor(hupsa$epi, levels = desired_order)
Idents(hupsa) <- hupsa$epi
DimPlot(hupsa, label = T)

DimPlot2(hupsa,reduction = "umap", group.by = "epi",theme = theme_umap_arrows())+
  scale_x_continuous(limits = c(-12, 15)) +
  scale_y_continuous(limits = c(-12, 15))

gene_list <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/4.1 Add_gene_signature/lineage_240222_aligned.rds")
for(gn in names(gene_list)){
  hupsa <- AddModuleScore(hupsa, 
                              features = list(gene_list[[gn]]), 
                              assay = "RNA",
                              name = gn)
}
DotPlot(object = hupsa, features = paste(names(gene_list), "1", sep=""), scale.by = "size")+
  scale_colour_gradient2(high = "#3F0A3F", mid = "#B6CDE3", low = "white", midpoint = 0)+ 
  RotatedAxis()


#-----------
data <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/epi_rough(241024).rds")
data <- UpdateSeuratObject(data)

### start
anchor <- FindTransferAnchors(
  reference = hupsa,
  query = data,
  reference.reduction = "pca",
  normalization.method = "LogNormalize",
  dims = 1:50
)

gc()

data=MapQuery(
  anchorset = anchor,
  query = data,
  reference = hupsa,
  reference.reduction = "pca",
  refdata = list(
    epi = "epi",
    cell_type = "cell_type",
    histo = "histo"
  ),
  reduction.model = "umap"
)

colors<-c(
  "AdPCa-1" = "#6E2A26",
  "NEPC-1" = "#A9542A",
  "AdPCa-3" = "#C9C36B",
  "Progenitor-1" = "#D1A0B0",
  "AdPCa-6" = "#42644A",
  "AdPCa-2" = "#6B997E",
  "AdPCa-7" = "#EFC195",
  "AdPCa-5" = "#5E8EC7",
  "NEPC-2" = "#E95793",
  "AdPCa-4" = "#9E81BC",
  "Progenitor-2"   = "#A55089"
)

DimPlot(data,reduction = "ref.umap", group.by = "predicted.epi", pt.size = 1,split.by = "Sample_ID")+
  coord_cartesian(xlim = c(-12, 15), ylim = c(-15, 15)) 

DimPlot(hupsa, group.by = "epi", cols = col)+
  coord_cartesian(xlim = c(-12, 15), ylim = c(-15, 15)) 

DimPlot(hupsa, group.by = "epi", label = T)+
  coord_cartesian(xlim = c(-12, 15), ylim = c(-15, 15)) 

DimPlot(data,reduction = "ref.umap", group.by = "predicted.cell_type")+
  coord_cartesian(xlim = c(-12, 15), ylim = c(-15, 15)) 



DimPlot(data,reduction = "ref.umap", group.by = "Sample_ID", pt.size = 1)+
  coord_cartesian(xlim = c(-15, 15), ylim = c(-15, 15)) 



colors<-c(
  "AdPCa-1" = "#6A332F",
  "AdPCa-2" = "#987C73",
  "AdPCa-3" = "#B36C4B",
  "AdPCa-5" = "#3C492E",
  "NEPC-1" = "#6E78B5",
  "NEPC-2" = "#423F52",
  "Progenitor-1" = "#573566"
)
desired_order <- c("sgNT", "sgNLN", "sgNLN_WNT")
data$Sample_ID <- factor(data$Sample_ID, levels = desired_order)
DimPlot2(data,reduction = "ref.umap", group.by = "predicted.epi", cols = colors, pt.size = 0.5, split.by = "Sample_ID",
         shuffle = T, theme = theme_umap_arrows())+
  scale_x_continuous(limits = c(-12, 15)) +
  scale_y_continuous(limits = c(-12, 15))


colors<-c(
  "AdPCa-1" = "#6E2A26",
  "NEPC-1" = "#A9542A",
  "AdPCa-3" = "#C9C36B",
  "Progenitor-1" = "#D1A0B0",
  "AdPCa-6" = "#42644A",
  "AdPCa-2" = "#6B997E",
  "AdPCa-7" = "#EFC195",
  "AdPCa-5" = "#5E8EC7",
  "NEPC-2" = "#bf812d",
  "AdPCa-4" = "#9E81BC",
  "Progenitor-2"   = "#A55089"
)

DimPlot2(hupsa,reduction = "umap", group.by = "epi",cols = colors,theme = theme_umap_arrows())+
  scale_x_continuous(limits = c(-12, 15)) +
  scale_y_continuous(limits = c(-12, 15))


VolcanoPlot(data, 
            ident.1 = "B cell",
            ident.2 = "CD8 T cell")

# saveRDS(data, "/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/epi_ref(250506).rds")

#############################################
meta <- FetchData(data, vars = c("Sample_ID", "predicted.epi"))
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


ggplot(meta_percentages, aes(x = Sample_ID, y = percentage, fill = cluster)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  labs(y = "Percentage (%)", x = "Sample_ID", title = "Cell Type Proportions by ID") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))+
  coord_flip()

