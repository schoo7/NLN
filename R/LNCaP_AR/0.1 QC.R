rm(list=ls())
gc()
# Load packages and functions -----------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
#remotes::install_github('https://github.com/ekernf01/DoubletFinder', force = T)
library(DoubletFinder)

#-------------------------------------------------------------------------------
# perform cleaning and get the mitochondrial exp percentage
qc <- function(sobj){
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  
  return(sobj)
}
#-------------------------------------------------------------------------------
# Plot FeatureScatter
plotfeat <- function(sobj){
  plot1 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme(legend.position = "none")
  plot2 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(legend.position = "none")
  plot1 + plot2 
}
#-------------------------------------------------------------------------------
# Standard Seurat workflow, needed for the DoubletFinder
stdwkfl <- function(sobj){
  sobj <- NormalizeData(sobj)
  sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
  sobj <- ScaleData(sobj)
  sobj <- RunPCA(sobj)
  sobj <- RunUMAP(sobj, dims = 1:20)
  sobj <- FindNeighbors(sobj, dims = 1:20)
  sobj <- FindClusters(sobj, resolution = 1)
  
  return(sobj)
}
#-------------------------------------------------------------------------------
rundbfinder <- function(sobj){
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.lst <- paramSweep_v3(sobj, PCs = 1:15, sct = FALSE)
  sweep <- summarizeSweep(sweep.lst, GT = FALSE)
  
  pdf(file="/work/SCCC/s228508/grey_Mu_lab/choushi/bcmvn/bcmvn_figure.pdf", width=10, height=6)
  bcmvn <- find.pK(sweep)
  dev.off()
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(sobj@meta.data$seurat_clusters)  # ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(sobj@meta.data))      # Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  # set pK to the pK value of max BCmetric
  pK = as.numeric(as.character(bcmvn[bcmvn$BCmetric == max(bcmvn$BCmetric),]$pK))
  sobj.1 = doubletFinder_v3(sobj, PCs = 1:15, pK = pK, nExp = nExp_poi)
  sobj.2 <- doubletFinder_v3(sobj, PCs = 1:15, pK = pK, nExp = nExp_poi.adj)
  
  ## Get the column names for each ----------------------------------------------------------------
  db1 = paste("DF.classifications_0.25_", pK, "_", nExp_poi, sep = "")
  db2 = paste("DF.classifications_0.25_", pK, "_", nExp_poi.adj, sep = "")
  
  ## Rename the columns ----------------------------------------------------------------
  names(sobj.1@meta.data)[names(sobj.1@meta.data) == db1] <- "DB1"
  names(sobj.2@meta.data)[names(sobj.2@meta.data) == db2] <- "DB2"
  
  ## High confidence if Doublet in both params ----------------------------------------------------------------
  sobj.1@meta.data$DB2 = sobj.2@meta.data$DB2
  meta = sobj.1@meta.data %>% mutate(DF = ifelse((DB1 == "Doublet") & DB2 == "Doublet", "Doublet - High Confidence", 
                                                 ifelse((DB2 == "Doublet" & DB1 == "Singlet") | (DB2 == "Singlet" & DB1 == "Doublet"), "Doublet - Low Confidence", "Singlet")))
  
  ## Put it back into original seurat obj ----------------------------------------------------------------
  sobj@meta.data$DF = meta$DF
  return(sobj)
}
#-------------------------------------------------------------------------------
savenew <- function(sobj,dir){
  saveRDS(sobj, paste("/work/SCCC/s228508/grey_Mu_lab/choushi/after_qc_all_genes/", dir, sep = ""))
}
#-------------------------------------------------------------------------------
plotdensity <- function(sobj){
  p <- ggplot(sobj@meta.data, aes(x=nFeature_RNA)) + 
    geom_density()
  return(p)
}
#-------------------------------------------------------------------------------
get2madcutsoff <-function(sobj){
  res = median(sobj@meta.data$nFeature_RNA) - 2* mad(sobj@meta.data$nFeature_RNA)
  if (res < 0){
    res = 0
  }
  return(res)
}

# nFeature_RNA > threshold(get2madcutsoff)
#-------------------------------------------------------------------------------
diagnosticplots <- function(sobj, vline = NULL){
  p1 <- plotfeat(sobj)
  
  if (!is.null(vline)){
    p2 <- plotdensity(sobj) + geom_vline(xintercept = vline)
  } else {
    p2 <- plotdensity(sobj)
  }
  
  print(get2madcutsoff(sobj))
  
  return(p1+p2)
  
}







#############################################################################################################
setwd("/work/SCCC/s228508/grey_Mu_lab/choushi/seurat_obj_all_genes/")

#############################################################################################################
sobj <- readRDS("sgNT.rds")
sobj <- qc(sobj)

# pdf(file="/work/SCCC/s228508/integration/QC_plots/10879_F_diagnose.pdf", width=10, height=6)
diagnosticplots(sobj)
# dev.off()

# pdf(file="/work/SCCC/s228508/integration/QC_plots/10879_F_violin.pdf", width=10, height=6)
VlnPlot(object = sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)+ #gene numbers, depth of seq, percentage of mitochondrial genes
  geom_hline(yintercept = c(10,15,20), linetype = "dashed", color = "red")
# dev.off()

# Apply cutoff
# sobj <- subset(sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15)
sobj <- subset(sobj, subset = nFeature_RNA > 200 & percent.mt < 20) #15 is ok

# Run DoubletFinder and subset
sobj <- stdwkfl(sobj)

# pdf(file="/work/SCCC/s228508/integration/QC_plots/10879_F_plot.pdf", width=10, height=6)
FeaturePlot(sobj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt")) + DimPlot(sobj)
# dev.off()

sobj <- rundbfinder(sobj)
table(sobj$DF)
# sobj <- SetIdent(sobj, value = sobj@meta.data$seurat_clusters)
# plot1 <- DimPlot(sobj)
# sobj <- SetIdent(sobj, value = sobj@meta.data$DF)
# plot2 <- DimPlot(sobj)
# plot1 + plot2
sobj <- subset(sobj, subset = DF == "Doublet - High Confidence", invert = T) #it keeps cells that do not match the condition

# Save the QC'd seurat
# sobj <- DietSeurat(sobj)
savenew(sobj,"sgNT.rds")


#------------------------------------
sobj <- readRDS("sgNLN.rds")
sobj <- qc(sobj)

# pdf(file="/work/SCCC/s228508/integration/QC_plots/10879_F_diagnose.pdf", width=10, height=6)
diagnosticplots(sobj)
# dev.off()

# pdf(file="/work/SCCC/s228508/integration/QC_plots/10879_F_violin.pdf", width=10, height=6)
VlnPlot(object = sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)+ 
  geom_hline(yintercept = c(10,15,20), linetype = "dashed", color = "red")
# dev.off()

# Apply cutoff
# sobj <- subset(sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15)
sobj <- subset(sobj, subset = nFeature_RNA > 200 & percent.mt < 20) #15 is ok

# Run DoubletFinder and subset
sobj <- stdwkfl(sobj)

# pdf(file="/work/SCCC/s228508/integration/QC_plots/10879_F_plot.pdf", width=10, height=6)
FeaturePlot(sobj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt")) + DimPlot(sobj)
# dev.off()

sobj <- rundbfinder(sobj)
table(sobj$DF)

sobj <- subset(sobj, subset = DF == "Doublet - High Confidence", invert = T) 

savenew(sobj,"sgNLN.rds")





#------------------------------------
sobj <- readRDS("sgNLN_KIF.rds")
sobj <- qc(sobj)

# pdf(file="/work/SCCC/s228508/integration/QC_plots/10879_F_diagnose.pdf", width=10, height=6)
diagnosticplots(sobj)
# dev.off()

# pdf(file="/work/SCCC/s228508/integration/QC_plots/10879_F_violin.pdf", width=10, height=6)
VlnPlot(object = sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)+ 
  geom_hline(yintercept = c(10,15,20), linetype = "dashed", color = "red")
# dev.off()

# Apply cutoff
# sobj <- subset(sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15)
sobj <- subset(sobj, subset = nFeature_RNA > 200 & percent.mt < 20) #15 is ok

# Run DoubletFinder and subset
sobj <- stdwkfl(sobj)

# pdf(file="/work/SCCC/s228508/integration/QC_plots/10879_F_plot.pdf", width=10, height=6)
FeaturePlot(sobj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt")) + DimPlot(sobj)
# dev.off()

sobj <- rundbfinder(sobj)
table(sobj$DF)

sobj <- subset(sobj, subset = DF == "Doublet - High Confidence", invert = T) 

savenew(sobj,"sgNLN_KIF.rds")



#------------------------------------
sobj <- readRDS("sgNLN_WNT.rds")
sobj <- qc(sobj)

# pdf(file="/work/SCCC/s228508/integration/QC_plots/10879_F_diagnose.pdf", width=10, height=6)
diagnosticplots(sobj)
# dev.off()

# pdf(file="/work/SCCC/s228508/integration/QC_plots/10879_F_violin.pdf", width=10, height=6)
VlnPlot(object = sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)+ 
  geom_hline(yintercept = c(10,15,20), linetype = "dashed", color = "red")
# dev.off()

# Apply cutoff
# sobj <- subset(sobj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 15)
sobj <- subset(sobj, subset = nFeature_RNA > 200 & percent.mt < 20) #15 is ok

# Run DoubletFinder and subset
sobj <- stdwkfl(sobj)

# pdf(file="/work/SCCC/s228508/integration/QC_plots/10879_F_plot.pdf", width=10, height=6)
FeaturePlot(sobj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt")) + DimPlot(sobj)
# dev.off()

sobj <- rundbfinder(sobj)
table(sobj$DF)

sobj <- subset(sobj, subset = DF == "Doublet - High Confidence", invert = T) 

savenew(sobj,"sgNLN_WNT.rds")





