rm(list = ls())

library(Seurat)
library(Matrix)
library(tidyverse)
library(org.Hs.eg.db)
library(biomaRt)

mart <- useMart("ensembl","hsapiens_gene_ensembl")

# listDatasets(useMart("ensembl"))$dataset
# tmp <- listAttributes(mart)
# DT::datatable(tmp)

anno <- getBM(values = gene_table$symbols,
              mart = mart,
              filters ="ensembl_gene_id",
              attributes = c("ensembl_gene_id", "external_gene_name"))

anno$external_gene_name[anno$external_gene_name == ""] <- NA
sum(is.na(anno))

anno_clean <- na.omit(anno, cols = "external_gene_name")
rownames(anno_clean) <- anno_clean$ensembl_gene_id

#------------------
count_matrix <- readMM("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNT/count_matrix.mtx")
all_genes <- read.csv("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNT/all_genes.csv", header = T)
count_matrix <- t(count_matrix)
head(count_matrix)
dim(count_matrix)

rownames(count_matrix) <- all_genes$gene_id
colnames(count_matrix) <- paste0("Cell_", 1:ncol(count_matrix))
# count_matrix <- data.frame(count_matrix)


shareID <- intersect(rownames(anno_clean), rownames(count_matrix))
count_matrix <- count_matrix[shareID,]

exprSet_aggregated <- aggregate(count_matrix, by = list(anno_clean$external_gene_name), FUN = max)
exprSet_aggregated <- column_to_rownames(exprSet_aggregated, 'Group.1')
dim(exprSet_aggregated)

seurat_obj <- CreateSeuratObject(counts = exprSet_aggregated, project = "sgNT", min.cells = 3,min.features = 100)

rownames_to_filter <- rownames(exprSet_aggregated)
rows_to_remove <- grepl("_", rownames_to_filter)
exprSet_filtered <- exprSet_aggregated[!rows_to_remove, ]
dim(exprSet_filtered)


#-----------------------------
# sgNT
count_matrix <- readMM("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNT/count_matrix.mtx")
count_matrix <- t(count_matrix)
head(count_matrix)
dim(count_matrix)

# rownames(count_matrix) <- paste0("Gene_", 1:nrow(count_matrix))
colnames(count_matrix) <- paste0("Cell_", 1:ncol(count_matrix))
count_matrix <- data.frame(count_matrix)

all_genes <- read.csv("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNT/all_genes.csv", header = T)
gene_table <- data.frame(cbind(genes=rownames(count_matrix), symbols = all_genes$gene_id))

# table(duplicated(gene_table$symbols))




exprSet_aggregated <- aggregate(count_matrix, by = list(gene_table$symbols), FUN = sum)
exprSet_aggregated <- column_to_rownames(exprSet_aggregated, 'Group.1')
dim(exprSet_aggregated)

rownames_to_filter <- rownames(exprSet_aggregated)
rows_to_remove <- grepl("^ENS", rownames_to_filter) | grepl("_", rownames_to_filter)
exprSet_filtered <- exprSet_aggregated[!rows_to_remove, ]
dim(exprSet_filtered)

# options(Seurat.object.assay.version = 'v5')

exprSet_sparse <- as(as.matrix(exprSet_filtered), "dgCMatrix")
seurat_obj <- CreateSeuratObject(counts = exprSet_sparse, project = "sgNT", min.cells = 3,min.features = 100)
saveRDS(seurat_obj, "/work/SCCC/s228508/grey_Xiaoling/choushi/seurat_obj/sgNT.rds")

#########################
# sgNLN
count_matrix <- readMM("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNLN/count_matrix.mtx")
count_matrix <- t(count_matrix)
head(count_matrix)
dim(count_matrix)

rownames(count_matrix) <- all_genes$gene_name
colnames(count_matrix) <- paste0("Cell_", 1:ncol(count_matrix))
count_matrix <- data.frame(count_matrix)

all_genes <- read.csv("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNLN/all_genes.csv", header = T)
gene_table <- data.frame(cbind(genes=rownames(count_matrix), symbols = all_genes$gene_name))

table(duplicated(gene_table$symbols))

exprSet_aggregated <- aggregate(count_matrix, by = list(gene_table$symbols), FUN = sum)
exprSet_aggregated <- column_to_rownames(exprSet_aggregated, 'Group.1')
dim(exprSet_aggregated)

rownames_to_filter <- rownames(exprSet_aggregated)
rows_to_remove <- grepl("^ENS", rownames_to_filter) | grepl("_", rownames_to_filter)
exprSet_filtered <- exprSet_aggregated[!rows_to_remove, ]
dim(exprSet_filtered)

# options(Seurat.object.assay.version = 'v5')

exprSet_sparse <- as(as.matrix(exprSet_filtered), "dgCMatrix")
seurat_obj <- CreateSeuratObject(counts = exprSet_sparse, project = "sgNLN", min.cells = 3,min.features = 100)
saveRDS(seurat_obj, "/work/SCCC/s228508/grey_Xiaoling/choushi/seurat_obj/sgNLN.rds")

#########################
# sgNLN_KIF
count_matrix <- readMM("/work/SCCC/s228508/grey_Xiaoling/choushi/raw_data/sgNLN_KIF/count_matrix.mtx")
count_matrix <- t(count_matrix)
head(count_matrix)
dim(count_matrix)

rownames(count_matrix) <- paste0("Gene_", 1:nrow(count_matrix))
colnames(count_matrix) <- paste0("Cell_", 1:ncol(count_matrix))
count_matrix <- data.frame(count_matrix)

all_genes <- read.csv("/work/SCCC/s228508/grey_Xiaoling/choushi/raw_data/sgNLN_KIF/all_genes.csv", header = T)
gene_table <- data.frame(cbind(genes=rownames(count_matrix), symbols = all_genes$gene_name))

table(duplicated(gene_table$symbols))

exprSet_aggregated <- aggregate(count_matrix, by = list(gene_table$symbols), FUN = sum)
exprSet_aggregated <- column_to_rownames(exprSet_aggregated, 'Group.1')
dim(exprSet_aggregated)

rownames_to_filter <- rownames(exprSet_aggregated)
rows_to_remove <- grepl("^ENS", rownames_to_filter) | grepl("_", rownames_to_filter)
exprSet_filtered <- exprSet_aggregated[!rows_to_remove, ]
dim(exprSet_filtered)

# options(Seurat.object.assay.version = 'v5')

exprSet_sparse <- as(as.matrix(exprSet_filtered), "dgCMatrix")
seurat_obj <- CreateSeuratObject(counts = exprSet_sparse, project = "sgNLN_KIF", min.cells = 3,min.features = 100)
saveRDS(seurat_obj, "/work/SCCC/s228508/grey_Xiaoling/choushi/seurat_obj/sgNLN_KIF.rds")

########################
# sgNLN_WNT
count_matrix <- readMM("/work/SCCC/s228508/grey_Xiaoling/choushi/raw_data/sgNLN_WNT/count_matrix.mtx")
count_matrix <- t(count_matrix)
head(count_matrix)
dim(count_matrix)

rownames(count_matrix) <- paste0("Gene_", 1:nrow(count_matrix))
colnames(count_matrix) <- paste0("Cell_", 1:ncol(count_matrix))
count_matrix <- data.frame(count_matrix)

all_genes <- read.csv("/work/SCCC/s228508/grey_Xiaoling/choushi/raw_data/sgNLN_WNT/all_genes.csv", header = T)
gene_table <- data.frame(cbind(genes=rownames(count_matrix), symbols = all_genes$gene_name))

table(duplicated(gene_table$symbols))

exprSet_aggregated <- aggregate(count_matrix, by = list(gene_table$symbols), FUN = sum)
exprSet_aggregated <- column_to_rownames(exprSet_aggregated, 'Group.1')
dim(exprSet_aggregated)

rownames_to_filter <- rownames(exprSet_aggregated)
rows_to_remove <- grepl("^ENS", rownames_to_filter) | grepl("_", rownames_to_filter)
exprSet_filtered <- exprSet_aggregated[!rows_to_remove, ]
dim(exprSet_filtered)

# options(Seurat.object.assay.version = 'v5')

exprSet_sparse <- as(as.matrix(exprSet_filtered), "dgCMatrix")
seurat_obj <- CreateSeuratObject(counts = exprSet_sparse, project = "sgNLN_WNT", min.cells = 3,min.features = 100)
saveRDS(seurat_obj, "/work/SCCC/s228508/grey_Xiaoling/choushi/seurat_obj/sgNLN_WNT.rds")


#------------
# Install SeuratDisk if not already installed
# remotes::install_github("mojaveazure/seurat-disk")
# 
# BiocManager::install("SeuratDisk")
# # Load the libraries
# library(SeuratDisk)
# library(Seurat)
# 
# # Convert the .h5ad file to .h5Seurat
# Convert("path_to_file/anndata.h5ad", dest = "h5seurat", overwrite = TRUE)
# # Load the .h5Seurat file into a Seurat object
# seurat_obj <- LoadH5Seurat("path_to_file/anndata.h5Seurat")
# # Check the Seurat object
# seurat_obj





