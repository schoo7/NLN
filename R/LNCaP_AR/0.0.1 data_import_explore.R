rm(list = ls())

seurat.primary <- readRDS("/work/SCCC/s228508/grey_Mu_lab/su/diet_seurat_avereps/sgIWS1_Veh.rds")

seurat.primary <- sobj

DefaultAssay(seurat.primary) <- "RNA"
seurat.primary <- NormalizeData(seurat.primary)
seurat.primary <- FindVariableFeatures(seurat.primary, nfeatures = 2000)

all.genes <- rownames(seurat.primary)
seurat.primary <- ScaleData(seurat.primary, features = all.genes, verbose = T)
seurat.primary <- RunPCA(seurat.primary,npcs = 50, verbose = T)
ElbowPlot(seurat.primary, ndims = 50)

set.seed(123)
pcSelect = 20
seurat.primary <- FindNeighbors(seurat.primary, reduction = "pca", dims = 1:pcSelect)
seurat.primary <- FindClusters(seurat.primary, resolution = 0.5) 
seurat.primary <- RunUMAP(seurat.primary, reduction = "pca", dims = 1:pcSelect)

# saveRDS(seurat.harmony, "/work/SCCC/s228508/grey_KR/primary_harmony_c16(240802).rds")

DimPlot(seurat.primary, reduction = "umap", label = T, raster= F) 






#-------------
library(SeuratDisk)
library(rhdf5)
# Convert the .h5ad file to .h5Seurat format
Convert("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNT/anndata.h5ad", dest = "h5seurat")

# Load the .h5Seurat file as a Seurat object
seurat_obj <- LoadH5Seurat("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNT/anndata.h5seurat", meta.data = FALSE, misc = FALSE)


h5ls("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNT/anndata.h5ad")



#-----------
library(zellkonverter)
library(SingleCellExperiment)

ad <- readH5AD("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNLN/anndata.h5ad")
assays(ad)  # should show the count matrices
dim(ad)  # Number of genes x cells

rownames(ad)
colnames(ad)

colData(ad)
rowData(ad)

raw_counts <- assay(ad, "counts")  # May vary, check what’s available
head(raw_counts)
# List available assays
assayNames(ad)

# Check available layers
ad@int_metadata$h5ad$layers

normalized_data <- assay(ad, "X")  # Assuming "X" contains normalized data
# Convert log-normalized data back to counts (assuming log1p transformation)
raw_counts_approx <- expm1(normalized_data)


adata_Seurat <- as.Seurat(ad, counts = "X", data = NULL)

# diet_obj <- CreateSeuratObject(counts = adata_Seurat@assays[["originalexp"]]@counts)

count_table <- adata_Seurat@assays[["originalexp"]]@counts
count_table <- data.frame(count_table)
count_table <- data.frame(count_table)
# a <- readRDS("/work/SCCC/s228508/grey_Mu_lab/choushi/diet_seurat/sgNT.rds")

count_table2 <- 2^count_table - 1
count_table2 <- data.frame(count_table2)
# If needed, replace the logcounts assay with raw counts
assay(sce, "counts") <- counts_matrix





#----------
library("Seurat")
library("anndata")
options(Seurat.object.assay.version = 'v3')

print("Convert from Scanpy to Seurat...")
data <- read_h5ad("/work/SCCC/s228508/grey_Mu_lab/choushi/raw_data/sgNT/anndata.h5ad")
data <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
print(str(data))

















