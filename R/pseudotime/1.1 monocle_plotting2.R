rm(list=ls())
library(monocle3)
library(ggplot2)
library(devtools)
library(dplyr)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

# setwd('/work/SCCC/s228508/10.0 pseudotime_trajectory/')
# cds2 = readRDS('monocle.Rdata')
# expr2 = as.matrix(exprs(cds2))
cds = readRDS('/gpfs/gibbs/project/mu/qx75/Grey/10.0 pseudotime_trajectory/monocle.Rdata')
expr = as.matrix(exprs(cds))

head(expr)[1:5,1:5]

primary.tumor.epi <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/1.1.2 tumor_epi_reintegration/tepi_harmomy_p30r1.5_grouping_epi_rough.rds")

meta <- primary.tumor.epi@meta.data
expr = as.matrix(primary.tumor.epi@assays[["RNA"]]@counts)

########## custom plots ##########
gene_list <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/4.1 Add_gene_signature/lineage_240222_aligned.rds")
gene_list_selected = gene_list[[1]]

cluster_colors = setNames(brewer.pal(name="Set1", n=4), sort(unique(meta$epi_rough)))
pseudotime_colors = colorRamp2(seq(0, 15, length = 9), rev(brewer.pal(name="YlGnBu", n=9)))

topanno = HeatmapAnnotation(
    Pseudotime = pseudotime(cds)[order(pseudotime(cds))],
    Cluster = meta$epi_rough[order(pseudotime(cds))],
    col = list(Pseudotime=pseudotime_colors, Cluster=cluster_colors))


genes = unique(gene_list_selected)
print(genes)  # To view the genes
# genes = c("IL1B","VIM","KLK2")
genes = genes[genes %in% row.names(expr)]

plot_data = expr[genes, order(pseudotime(cds))] # Orders the expression data by pseudotime.
plot_data = t(apply(plot_data,1,function(x){smooth.spline(x,df=3)$y})) #Smooths the expression values for each gene across pseudotime using a smoothing spline
plot_data = t(apply(plot_data,1,function(x){(x-mean(x))/sd(x)})) #Standardizes (z-scores) the smoothed expression values to make gene expression comparable across genes
# psut = pseudotime(cds)/max(pseudotime(cds)) # Normalizes the pseudotime values to a [0, 1] range

# pdf(paste('pseudotime_', gl, '_top.pdf',sep=''), height=length(genes)/4)
g = Heatmap(
      plot_data,
      col = colorRamp2(seq(from=-2,to=2,length=9), rev(brewer.pal(11, "Spectral"))),
      show_row_names = TRUE,
      show_column_names  = FALSE,
      row_title_rot = 0,
      top_annotation  = topanno,
      cluster_row_slices = FALSE,
      cluster_columns = FALSE,
      cluster_row = T,
      use_raster =F)

