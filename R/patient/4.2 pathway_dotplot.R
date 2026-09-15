rm(list=ls())
library(GSEABase)
library(msigdbr)
library(RColorBrewer)
library(Seurat)
library(data.table)
library(colorspace)  
library(ggplot2)
library(dplyr)
library(harmony)
library(cowplot)
library(ComplexHeatmap)
library(circlize)

tumor.epi <- readRDS("/work/SCCC/s228508/choushi/tumor.epi_NLN_class(240511).rds")
colors <- c("NLN_low"="#D6606B", "NLN_high"="#4575b4")
DimPlot(tumor.epi, label = T, raster = F,cols = colors) 

DefaultAssay(tumor.epi) <- "RNA"

#--------------------
endocytosis <- readRDS("/work/SCCC/s228508/choushi/Endocytosis.rds")

#load HAKKMARK
m_df<- msigdbr(species = "human",  category = "H" )
geneSets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

sig <- readRDS("/work/SCCC/s228508/choushi/Dotplot_sig1.rds")

for(gn in names(sig)){
  tumor.epi <- AddModuleScore(tumor.epi, 
                              features = list(sig[[gn]]), 
                              assay = "RNA",
                              name = gn)
}

# meta <- FetchData(tumor.epi, vars = c("GOBP_CLATHRIN_DEPENDENT_ENDOCYTOSIS1", "HALLMARK_WNT_BETA_CATENIN_SIGNALING1","NLN_class"))
# colnames(meta) <- c("GOBP_CLATHRIN_DEPENDENT_ENDOCYTOSIS1", "HALLMARK_WNT_BETA_CATENIN_SIGNALING1","NLN_class")
# wilcox.test(HALLMARK_WNT_BETA_CATENIN_SIGNALING1 ~ NLN_class, data = meta, paired=FALSE)

pdf("/work/SCCC/s228508/choushi/figs/sig1.pdf", width=12, height=8)
DotPlot(object = tumor.epi, features = paste(names(sig), "1", sep=""), cols = c("#FFF1ED", "#B00015"), scale.by = "size")+
  scale_size(range = c(0.5,6))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()



