
# Example
reticulate::py_install("palantir", method = "pip", pip = TRUE)

# Download the example Seurat Object with myeloid cells
mye_small <- readRDS(url("https://zenodo.org/records/10944066/files/pbmc10k_mye_small_velocyto.rds", "rb"))

# Compute diffusion map
mye_small <- Palantir.RunDM(mye_small)

# Visualize the first two diffusion map dimensions
DimPlot2(mye_small, reduction = "ms")


# Calculate pseudotime with a specified start cell
mye_small <- Palantir.Pseudotime(mye_small, start_cell = "sample1_GAGAGGTAGCAGTACG-1")

# Store pseudotime results in meta.data for easy plotting
ps <- mye_small@misc$Palantir$Pseudotime
colnames(ps)[3:4] <- c("fate1", "fate2")
mye_small@meta.data[,colnames(ps)] <- ps

# Visualize pseudotime and cell fates
DimPlot2(
  mye_small,
  features = colnames(ps),
  reduction = "ms",
  cols = list(Entropy = "D"),
  theme = NoAxes())


###########################################
tumor.epi <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/epi_rough(241024).rds")

tumor.epi <- Palantir.RunDM(tumor.epi)

# Visualize the first two diffusion map dimensions
DimPlot2(tumor.epi, reduction = "ms")


# Calculate pseudotime with a specified start cell
mye_small <- Palantir.Pseudotime(tumor.epi, start_cell = "sample1_GAGAGGTAGCAGTACG-1")

# Store pseudotime results in meta.data for easy plotting
ps <- mye_small@misc$Palantir$Pseudotime
colnames(ps)[3:4] <- c("fate1", "fate2")
mye_small@meta.data[,colnames(ps)] <- ps

# Visualize pseudotime and cell fates
DimPlot2(
  mye_small,
  features = colnames(ps),
  reduction = "ms",
  cols = list(Entropy = "D"),
  theme = NoAxes())

