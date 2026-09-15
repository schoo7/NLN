rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(GSEABase)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(enrichplot)
# remotes::install_github('https://github.com/NicolasH2/gggsea', force = T)
library(gggsea)
library(ggsci)

tumor.epi <- readRDS("/work/SCCC/s228508/choushi/tumor.epi_NLN_class(240511).rds")
colors <- c("NLN_low"="#D6606B", "NLN_high"="#4575b4")
DimPlot(tumor.epi, label = T, raster = F,cols = colors) 

# deg <- FindMarkers(tumor.epi,ident.1 = "NLN_low",ident.2 = "NLN_high",min.pct = 0.25,logfc.threshold = 0.01)
# write.csv(deg, "/work/SCCC/s228508/choushi/deg2.csv")

# deg <- read.table("/work/SCCC/s228508/choushi/deg.csv", sep = ",", header = T)
deg <- read.table("/work/SCCC/s228508/choushi/deg2.csv", sep = ",", header = T)
id = deg$avg_log2FC
names(id) = deg$X
id = sort(id,decreasing = T)
id

# load gene sets
#--------
GO_gmt <- read.gmt("/work/SCCC/s228508/5.2.1 CD4/GO_genesets.gmt")
GOALL_csv <- read.table("/work/SCCC/s228508/5.2.1 CD4/GOALL_genesets.csv", sep = ",")

endocytosis <- readRDS("/work/SCCC/s228508/choushi/Endocytosis.rds")
endocytosis_gmt <- read.gmt("/work/SCCC/s228508/choushi/endocytosis.gmt")
endocytosis_selected_gmt <- read.gmt("/work/SCCC/s228508/choushi/endocytosis_selected.gmt")
endocytosis_selected_gmt <- subset(endocytosis_selected_gmt, endocytosis_selected_gmt$term %in% 
                                     c("GOBP_VESICLE_TRANSPORT_ALONG_MICROTUBULE"
                                       #"GOBP_EXOCYTIC_PROCESS",
                                       #"GOBP_POSITIVE_REGULATION_OF_EXOSOMAL_SECRETION",
                                       #"GOBP_CLATHRIN_COATED_VESICLE_CARGO_LOADING"
                                     ))

hallmark_gmt <- read.gmt("/work/SCCC/s228508/5.2.1 CD4/HALLMARK.gmt")

hallmark_gmt <- subset(hallmark_gmt, hallmark_gmt$term %in% c("HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_WNT_BETA_CATENIN_SIGNALING"))

# run GSEA
egmt <- clusterProfiler::GSEA(id, TERM2GENE = endocytosis_gmt,
                              pvalueCutoff = 1, minGSSize = 1, maxGSSize = 1000, eps = 0, nPermSimple = 10000, pAdjustMethod = "fdr")

gsea_result <- egmt@result
gsea_result2 <- subset(gsea_result, p.adjust < 0.05)
gsea_result2 <- gsea_result2[order(gsea_result2$NES, decreasing = T),]

write.table(gsea_result2, "/work/SCCC/s228508/choushi/figs/gsea_GO(NLNLOWvsNLNhigh).csv", row.name = F, sep =",")

# # run fgsea
# fgseaRes<- fgseaMultilevel(pathways = endocytosis, stats = id, minSize = 1, maxSize = 1000)
# 
# g2 <- fgseaRes[fgseaRes$padj < 0.05,]
# g2 <- g2[order(g2$NES, decreasing = T),]
# 
# # check the top2 
# gseaplot2(egmt,geneSetID = c(1,2), pvalue_table = T, base_size = 5)



# visualization
#-----------------
# No.1
plotGseaTable(hallmark[g2$pathway],
              id,
              fgseaRes,
              gseaParam=0.5,
              colwidths = c(0.5, 0.2, 0.1, 0.1, 0.1))

# No.2
ggplot(fgseaRes %>% filter(padj < 0.05) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways") +
  theme_minimal() 

# No.3
color <- pal_simpsons()(16)

num1=1

pdf("/work/SCCC/s228508/choushi/figs/GSEA.pdf", width=5, height=4)
gseaplot2(egmt, geneSetID = rownames(gsea_result2)[1:num1],
          title = "GOBP_VESICLE_TRANSPORT_ALONG_MICROTUBULE",
          color = color[1:num1],
          base_size = 10,
          rel_heights = c(0.8,0.2,0.2),
          subplots = 1:3,
          pvalue_table = F,
          ES_geom = "line" #line or dot
)
dev.off()




