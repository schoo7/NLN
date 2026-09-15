rm(list=ls())

# library("BiocManager") 
# install dependencies
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'HDF5Array',
#                        'terra', 'ggrastr'))

# BiocManager::install('batchelor',force = TRUE)
# BiocManager::install('terra')

# remotes::install_github("rspatial/terra")

# devtools::install_github('cole-trapnell-lab/monocle3')
# devtools::install_github('junjunlab/ClusterGVis')

library(monocle3)
library(Seurat)
library(tidyverse)
library(ggplot2)
# library(ClusterGVis)
library(pheatmap)
library(reshape2)
library(viridis)

options(Seurat.object.assay.version = 'v3')

primary.tumor.epi <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/1.1.2 tumor_epi_reintegration/tepi_harmomy_p30r1.5_grouping_epi_rough.rds")
# DefaultAssay(primary.tumor.epi) <- "RNA"
# DimPlot(primary.tumor.epi)

data <- GetAssayData(primary.tumor.epi, assay = "RNA", slot = "counts")
cell_metadata <- primary.tumor.epi@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = row.names(data))
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

#----
# normalization
# cds <- preprocess_cds(cds)
cds <- preprocess_cds(cds, norm_method = c("log"), num_dim = 100) 
monocle3::plot_pc_variance_explained(cds)

#---
# batch effect
cds <- align_cds(cds, aligment_group = "Patient_ID")

#----
# umap
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA") # PCA is default
plot_cells(cds)

colnames(colData(cds))
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "epi_rough") + ggtitle("cds.umap") # epi_rough

#-----
# use the previous UMAP Embeddings
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(primary.tumor.epi, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

plot_cells(cds, reduction_method="UMAP", color_cells_by="epi_rough") + ggtitle('int.umap')

#-------
# clustering
cds <- cluster_cells(cds,cluster_method = "leiden", k = 40)#leiden,louvain
# cds <- cluster_cells(cds,cluster_method = "leiden", resolution=1e-5)#leiden,louvain
plot_cells(cds, color_cells_by = "cluster",
           show_trajectory_graph = F)

# plot_cells(cds, color_cells_by = "partition")

#-------
# trajectory inference
# With an argument use_partition we can specify if we want to learn a disjoint graph (value TRUE-default) in each partition
cds = learn_graph(cds, close_loop = T, verbose = T, use_partition = T)#,learn_graph_control = list(ncenter=1000))
# head(colData(cds))

p <- plot_cells(cds, color_cells_by = "epi_rough", 
               label_groups_by_cluster = F,               
               label_leaves = T, 
               label_branch_points = T,
               graph_label_size = 4,
               cell_size = 0.5,
               trajectory_graph_color = "grey28",
               trajectory_graph_segment_size = 0.75)
p

# ggsave("Trajectory.pdf", plot = p, width = 8, height = 6)




#######################
# define the root
# method #1
cds <- order_cells(cds, reduction_method = "UMAP")

plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_cell_groups = F,
           label_leaves = F, 
           label_branch_points = F,
           trajectory_graph_color = "white",
           show_trajectory_graph =F)


# saveRDS(cds, file = '/work/SCCC/s228508/10.0 pseudotime_trajectory/monocle(241020).Rdata')

#-------------
# method #2
# a helper function to identify the root principal points:
# set the specific cell type as the root in time_bin
get_earliest_principal_node  <- function(cds, time_bin="AR high"){
  
  cell_ids <- which(colData(cds)[, "epi_rough"] == time_bin)    
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex  
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])  
  root_pr_nodes <-    
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]    
  root_pr_nodes
}

cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds, color_cells_by = "pseudotime", 
           label_cell_groups = F,
           label_leaves = F, 
           label_branch_points = F)






#########################
# Find genes that change as a function of pseudotime
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores= 70)

# deg <- Track_genes %>% arrange(q_value) %>% filter(status == "OK") %>% head()

deg <- Track_genes %>% filter(q_value < 0.05) %>% arrange(-morans_I)
deg
# Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(deg, "/work/SCCC/s228508/10.0 pseudotime_trajectory/trajectory_genes(241020).csv", row.names = F)

Track_genes_sig <- deg %>%
  top_n(n=12, morans_I) %>%
  pull(gene_short_name) %>%
  as.character()

Track_genes_sig

levels(Idents(primary.tumor.epi))
lineage_cds <- cds[rowData(cds)$gene_short_name %in% Track_genes_sig,
                   colData(cds)$epi_rough %in% c("lineage high","AR low","AR high","basal")]
#lineage_cds <- order_cells(lineage_cds)

# fig1
plot_genes_in_pseudotime(cds["IL1B",],
                         color_cells_by = "pseudotime", # epi_rough
                         # min_expr = 0.1, 
                         ncol = 1)
# genes change in pseudotime

my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("KLK2","KLK3","TM4SF1","IL1B"))) 
cds_subset <- cds[my_genes,]
colData(cds_subset)


plot_genes_in_pseudotime(cds_subset, color_cells_by = "grouping")
plot_genes_in_pseudotime(cds_subset, color_cells_by = "epi_rough" )



# FeaturePlot
plot_cells(cds, genes= Track_genes_sig,
           cell_size=1, 
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)



# produce violin plots
plot_genes_violin(cds_subset, group_cells_by="cell_type", ncol=2)




###########################
# Extract the pseudotime values and add to the Seurat object
primary.tumor.epi <- AddMetaData(
  object = primary.tumor.epi,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

FeaturePlot(primary.tumor.epi, features = Track_genes_sig, pt.size = 0.1, order = T) & scale_color_viridis_c()





################################
# Extract pseudotime information
pseudotime <- pseudotime(cds)  
head(pseudotime)

genes <- c("KLK2","KLK3","TM4SF1","IL1B")

# Extract normalized gene expression matrix
raw_counts <- counts(cds)
Size_Factor <- colData(cds)$Size_Factor
normalized_data <- t(t(raw_counts) / Size_Factor)
log_normalized_data <- log1p(normalized_data)
head(log_normalized_data[, 1:5])

log_normalized_data <- data.frame(log_normalized_data)

expr_data <- data.frame(t(log_normalized_data[genes, ]))
expr_data$pseudotime <- pseudotime
expr_data$grouping <- primary.tumor.epi$grouping
expr_data$epi_rough <- primary.tumor.epi$epi_rough

expr_data_melt <- melt(expr_data, id.vars = c("pseudotime","grouping","epi_rough"), variable.name = "gene", value.name = "expression")

# Plot the data using ggplot2
ggplot(expr_data_melt, aes(x = pseudotime, y = expression)) +
  geom_smooth(aes(color = gene), method = "loess", span = 0.2) +
  facet_wrap(~ gene, scales = "free_y", nrow = 3) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Pseudotime", y = "Expression (min-to-max)")

# Plot the data using ggplot2, with different colors for different classes
ggplot(expr_data_melt, aes(x = pseudotime, y = expression, color = grouping)) +  # Correct the mapping for color here
  geom_smooth(method = "loess", span = 0.2) +
  facet_wrap(~ gene, scales = "free_y", nrow = 3) +
  theme_minimal() +
  labs(x = "Pseudotime", y = "Expression (min-to-max)", color = "grouping") +
  theme(legend.position = "right")


###########################
# extract pseudotime to seuratobject
primary.tumor.epi$pseudotime <- pseudotime(cds)
# FeaturePlot(primary.tumor.epi, features = "pseudotime")

normalized <- data.frame(GetAssayData(primary.tumor.epi, assay = "RNA", slot = "counts"))

expr_data <- data.frame(t(normalized[genes, ]))
# primary.tumor.epi_pseudotime <- data.frame(primary.tumor.epi$pseudotime)

# primary.tumor.epi_pseudotime_reordered <- primary.tumor.epi_pseudotime[order(rownames(primary.tumor.epi_pseudotime), 
#                                                                              rownames(expr_data)), ]

expr_data$pseudotime <- primary.tumor.epi$pseudotime 

# Melt the data for easier plotting
expr_data_melt <- melt(expr_data, id.vars = "pseudotime", variable.name = "gene", value.name = "expression")

# Plot the data using ggplot2
ggplot(expr_data_melt, aes(x = pseudotime, y = expression)) +
  geom_smooth(aes(color = gene), method = "loess", span = 0.2) +
  facet_wrap(~ gene, scales = "free_y", nrow = 3) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Pseudotime", y = "Expression (min-to-max)")



########################
# cell proportion
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, epi_rough, fill = epi_rough)) + geom_violin()

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(epi_rough, monocle3_pseudotime), fill = epi_rough)) + 
  geom_violin() + 
  scale_fill_manual(values = c("#FF6347", "#4682B4", "#3CB371", "#FFD700"))  # Customize your colors here

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(epi_rough, monocle3_pseudotime), fill = epi_rough)) + 
  geom_violin() + 
  scale_fill_viridis(discrete = TRUE)  # Use TRUE for categorical variables

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(epi_rough, monocle3_pseudotime), fill = epi_rough)) + 
  geom_violin() + 
  scale_fill_brewer(palette = "Set1")












