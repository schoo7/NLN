rm(list=ls())

library(Polychrome)
library(ggbeeswarm)
library(ggthemes)
library(slingshot)
library(SingleCellExperiment)
library(tidyverse)
library(RColorBrewer)

primary.tumor.epi <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/epi_rough(241024).rds")

sce <- as.SingleCellExperiment(primary.tumor.epi, assay = "RNA")
sce

reducedDims(sce)$Harmony <- Embeddings(primary.tumor.epi, reduction = "harmony")

sce <- slingshot(sce, 
                 clusterLabels = 'epi_rough',  
                 reducedDim = 'Harmony',  
                 start.clus= "c0",  
                 end.clus = NULL)
colnames(colData(sce))
head(colData(sce))

# this define the cluster color. You can change it with different color scheme.
my_color <- createPalette(length(levels(sce$seurat_clusters)), c("#010101", "#ff0000"), M=1000)
names(my_color) <- c("c0","c1","c2","c3")

my_color <- c("#d53e4f", "#bf812d", "#5ab4ac", "#4575b4")
names(my_color) <- c("0","1","2","3")
# Extract numeric or character columns only
slingshot_df <- as.data.frame(colData(sce)[, sapply(colData(sce), is.atomic)])

# re-order y-axis for better figure: This should be tailored with your own cluster names
# slingshot_df$ident = factor(slingshot_df$ident, levels=c(4,2,1,0,3,5,6))

ggplot(slingshot_df, aes(x = slingPseudotime_2, y = seurat_clusters, 
                         colour = seurat_clusters)) +
  geom_quasirandom(groupOnX = FALSE, size = 0.5) +  # adjust size here
  theme_classic() +
  xlab("First Slingshot pseudotime") +
  ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime") +
  scale_colour_manual(values = my_color)


#-------
plot(reducedDims(sce)$Harmony, 
     col = my_color[as.character(sce$seurat_clusters)], 
     pch = 16, 
     asp = 1,
     cex = 0.5)  # <-- adjust this value for smaller/larger dots

legend("bottomleft",legend = names(my_color[levels(sce$seurat_clusters)]),  
       fill = my_color[levels(sce$seurat_clusters)])

sce$seurat_clusters

lnes <- getLineages(reducedDim(sce,"Harmony"),
                    sce$seurat_clusters, start.clus = "0")
lnes@metadata$lineages

lines(SlingshotDataSet(lnes), lwd=2, type = 'lineages', col = c("black"))

lines(SlingshotDataSet(sce), lwd=2, col=brewer.pal(9,"Set1"))


#-------------------
canonical_WNT <- c("HDAC11","PPARD","KREMEN2","PRKCZ","CUL3","RHOA","KLHL12","FLNA","NOTCH1","NCSTN","HDAC5","FZD4","GNAI1","PSEN2","LRP5",
                   "ATP6AP2","WNT5A","WNT3A","CTNNB1","PORCN","FZD1","FZD5","JAG1","JAG2","LEF1","TCF7","LGR4","TCF7L2","ATOH1","KLF4")




#################
library(dplyr)

# Assuming your data is stored in a data frame called 'df'
slingshot_df <- slingshot_df %>%
  mutate(sling = coalesce(slingPseudotime_1, slingPseudotime_2))

tumor.epi <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/epi_rough(250429).rds")
#inferno,magma,plasma,viridis,cividis
tumor.epi$sling <- slingshot_df$sling

FeaturePlot(object = tumor.epi, reduction = "harmony", features = "sling", pt.size = 0.1, label = F, label.color = "white", label.size = 5, order = T)+ 
  scale_color_viridis_c(option = "inferno")

# scale_color_gradientn(colors = brewer.pal(n = 9, name = "YlGnBu"))
# scale_color_viridis_c()
# scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


colors <- c("c0"="#d53e4f", "c1"="#bf812d", "c2"="#5ab4ac", "c3"="#4575b4")
DimPlot2(tumor.epi, reduction = "harmony",theme = theme_umap_arrows(), label = F, cols = colors)

DimPlot(object = tumor.epi, reduction = "harmony", group.by = "epi_rough", cols = colors)


