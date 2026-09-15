rm(list = ls())

library(Seurat)
library(Matrix)
library(tidyverse)
library(org.Hs.eg.db)
library(biomaRt)
library(limma)

options(Seurat.object.assay.version = 'v3')

mart <- useMart("ensembl","hsapiens_gene_ensembl")

# listDatasets(useMart("ensembl"))$dataset
# tmp <- listAttributes(mart)
# DT::datatable(tmp)


##############################
# sgNT
count_matrix <- readMM("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNT/count_matrix.mtx")
all_genes <- read.csv("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNT/all_genes.csv", header = T)
count_matrix <- t(count_matrix)
head(count_matrix)
dim(count_matrix)

rownames(count_matrix) <- all_genes$gene_id
colnames(count_matrix) <- paste0("Cell_", 1:ncol(count_matrix))
# count_matrix <- data.frame(count_matrix)

#----
anno <- getBM(values = all_genes$gene_id,
              mart = mart,
              filters ="ensembl_gene_id",
              attributes = c("ensembl_gene_id", "external_gene_name"))

anno$external_gene_name[anno$external_gene_name == ""] <- NA
sum(is.na(anno))

anno_clean <- na.omit(anno, cols = "external_gene_name")
rownames(anno_clean) <- anno_clean$ensembl_gene_id

#----
shareID <- intersect(rownames(anno_clean), rownames(count_matrix))
count_matrix <- count_matrix[shareID,]

rownames(count_matrix) <- anno_clean$external_gene_name
count_matrix <- avereps(count_matrix)
dim(count_matrix)

seurat_obj <- CreateSeuratObject(counts = count_matrix, project = "sgNT", min.features = 100)
saveRDS(seurat_obj, "/work/SCCC/s228508/grey_Mu_lab/choushi/seurat_obj_all_genes/sgNT.rds")




##################################
# sgNLN
count_matrix <- readMM("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNLN/count_matrix.mtx")
all_genes <- read.csv("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNLN/all_genes.csv", header = T)
count_matrix <- t(count_matrix)
head(count_matrix)
dim(count_matrix)

rownames(count_matrix) <- all_genes$gene_id
colnames(count_matrix) <- paste0("Cell_", 1:ncol(count_matrix))
# count_matrix <- data.frame(count_matrix)

#----
anno <- getBM(values = all_genes$gene_id,
              mart = mart,
              filters ="ensembl_gene_id",
              attributes = c("ensembl_gene_id", "external_gene_name"))

anno$external_gene_name[anno$external_gene_name == ""] <- NA
sum(is.na(anno))

anno_clean <- na.omit(anno, cols = "external_gene_name")
rownames(anno_clean) <- anno_clean$ensembl_gene_id

#----
shareID <- intersect(rownames(anno_clean), rownames(count_matrix))
count_matrix <- count_matrix[shareID,]

rownames(count_matrix) <- anno_clean$external_gene_name
count_matrix <- avereps(count_matrix)
dim(count_matrix)

seurat_obj <- CreateSeuratObject(counts = count_matrix, project = "sgNLN", min.features = 100)
saveRDS(seurat_obj, "/work/SCCC/s228508/grey_Mu_lab/choushi/seurat_obj_all_genes/sgNLN.rds")



#########################
# sgNLN_KIF
count_matrix <- readMM("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNLN_KIF/count_matrix.mtx")
all_genes <- read.csv("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNLN_KIF/all_genes.csv", header = T)
count_matrix <- t(count_matrix)
head(count_matrix)
dim(count_matrix)

rownames(count_matrix) <- all_genes$gene_id
colnames(count_matrix) <- paste0("Cell_", 1:ncol(count_matrix))
# count_matrix <- data.frame(count_matrix)

#----
anno <- getBM(values = all_genes$gene_id,
              mart = mart,
              filters ="ensembl_gene_id",
              attributes = c("ensembl_gene_id", "external_gene_name"))

anno$external_gene_name[anno$external_gene_name == ""] <- NA
sum(is.na(anno))

anno_clean <- na.omit(anno, cols = "external_gene_name")
rownames(anno_clean) <- anno_clean$ensembl_gene_id

#----
shareID <- intersect(rownames(anno_clean), rownames(count_matrix))
count_matrix <- count_matrix[shareID,]

rownames(count_matrix) <- anno_clean$external_gene_name
count_matrix <- avereps(count_matrix)
dim(count_matrix)

seurat_obj <- CreateSeuratObject(counts = count_matrix, project = "sgNLN_KIF", min.features = 100)
saveRDS(seurat_obj, "/work/SCCC/s228508/grey_Mu_lab/choushi/seurat_obj_all_genes/sgNLN_KIF.rds")






#########################
# sgNLN_WNT
count_matrix <- readMM("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNLN_WNT/count_matrix.mtx")
all_genes <- read.csv("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNLN_WNT/all_genes.csv", header = T)
count_matrix <- t(count_matrix)
head(count_matrix)
dim(count_matrix)

rownames(count_matrix) <- all_genes$gene_id
colnames(count_matrix) <- paste0("Cell_", 1:ncol(count_matrix))
# count_matrix <- data.frame(count_matrix)

#----
anno <- getBM(values = all_genes$gene_id,
              mart = mart,
              filters ="ensembl_gene_id",
              attributes = c("ensembl_gene_id", "external_gene_name"))

anno$external_gene_name[anno$external_gene_name == ""] <- NA
sum(is.na(anno))

anno_clean <- na.omit(anno, cols = "external_gene_name")
rownames(anno_clean) <- anno_clean$ensembl_gene_id

#----
shareID <- intersect(rownames(anno_clean), rownames(count_matrix))
count_matrix <- count_matrix[shareID,]

rownames(count_matrix) <- anno_clean$external_gene_name
count_matrix <- avereps(count_matrix)
dim(count_matrix)

seurat_obj <- CreateSeuratObject(counts = count_matrix, project = "sgNLN_WNT", min.features = 100)
saveRDS(seurat_obj, "/work/SCCC/s228508/grey_Mu_lab/choushi/seurat_obj_all_genes/sgNLN_WNT.rds")



