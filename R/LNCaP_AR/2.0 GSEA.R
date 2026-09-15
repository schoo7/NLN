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

tumor.epi <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/epi_rough(250429).rds")

DefaultAssay(tumor.epi) <- "RNA" 
Idents(tumor.epi) <- "epi_rough"

markers <- FindMarkers(tumor.epi, ident.1 = "sgNLN", ident.2 = "sgNT", 
                       group.by = 'Sample_ID',
                       only.pos = F,
                       min.pct = 0.1,
                       logfc.threshold = 0.01)

id = markers$avg_log2FC
names(id) = rownames(markers)
id = sort(id,decreasing = T)
id

markers <- FindMarkers(tumor.epi, ident.1 = "c1", 
                       group.by = 'epi_rough',
                       only.pos = F,
                       min.pct = 0.1,
                       logfc.threshold = 0.01)

id = markers$avg_log2FC
names(id) = rownames(markers)
id = sort(id,decreasing = T)
id


# load gene sets
#--------
GOBP_gmt <- read.gmt("/work/SCCC/s228508/5.2.1 CD4/GOBP_genesets.gmt")
KEGG_gmt <- read.gmt("/work/SCCC/s228508/5.2.1 CD4/KEGG_genesets.gmt")
# CP_gmt <- read.gmt("/work/SCCC/s228508/5.2.1 CD4/CP_genesets.gmt")
# GO_gmt <- read.gmt("/work/SCCC/s228508/5.2.1 CD4/GO_genesets.gmt")
GOALL_csv <- read.table("/work/SCCC/s228508/5.2.1 CD4/GOALL_genesets.csv", sep = ",")
GOALL_csv2 <- read.table("/work/SCCC/s228508/7.0 Myeloid/all_GO_genesets2.csv", sep = ",")
# run GSEA
# GOBP
egmt <- clusterProfiler::GSEA(id, TERM2GENE = GOBP_gmt,
                              pvalueCutoff = 1, minGSSize = 10, maxGSSize = 1000, eps = 0, pAdjustMethod = "fdr")
# KEGG
set.seed(123)
egmt <- clusterProfiler::GSEA(id, TERM2GENE = KEGG_gmt,
                              pvalueCutoff = 1, minGSSize = 10, maxGSSize = 1000, eps = 0, pAdjustMethod = "fdr")
# GOALL2
set.seed(123)
egmt <- clusterProfiler::GSEA(id, TERM2GENE = GOALL_csv2,
                              pvalueCutoff = 1, minGSSize = 1, maxGSSize = 1000, eps = 0, nPermSimple = 1000, pAdjustMethod = "fdr") #nPermSimple = 10000, 
# GOALL
set.seed(123)
egmt <- clusterProfiler::GSEA(id, TERM2GENE = GOALL_csv,
                              pvalueCutoff = 1, minGSSize = 10, maxGSSize = 1000, eps = 0, nPermSimple = 1000, pAdjustMethod = "fdr") #nPermSimple = 10000, 

gsea_result <- egmt@result
gsea_result2 <- subset(gsea_result, p.adjust < 0.25)
gsea_result2 <- gsea_result2[order(gsea_result2$NES, decreasing = T),]

# write.table(gsea_result2, "/work/SCCC/s228508/7.0 Myeloid/GSEA/gsea_TAM_c06_inflammtory.csv", sep =",")

#----------------------------------
# run fgsea
GOALL_rds <- readRDS("/work/SCCC/s228508/7.0 Myeloid/GOALL_genesets(240814).rds")
ann_markers <- read.table("/work/SCCC/s228508/7.0 Myeloid/myeloid_rough_markers(240814).csv", sep = ",")
Endocytosis <- readRDS("/work/SCCC/s228508/choushi/Endocytosis.rds")

HALLMARK <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/7.0 Myeloid/HALLMARK(240816).rds")
for (i in names(HALLMARK)) {
  HALLMARK[[i]] <- unique(HALLMARK[[i]])
}

kegg <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/7.0 Myeloid/KEGG(240929).rds")
for (i in names(kegg)) {
  kegg[[i]] <- unique(kegg[[i]])
}

reactomes <- readRDS("/work/SCCC/s228508/7.0 Myeloid/reactomes(240914).rds")
for (i in names(reactomes)) {
  reactomes[[i]] <- unique(reactomes[[i]])
}

GOALL <- readRDS("/work/SCCC/s228508/7.0 Myeloid/GOALL_genesets(240814).rds")
for (i in names(GOALL)) {
  GOALL[[i]] <- unique(GOALL[[i]])
}

combined_list <- c(HALLMARK, kegg)  

fgseaRes<- fgseaMultilevel(pathways = combined_list, stats = id, minSize = 10, maxSize = 1000, nPermSimple = 1000)

g2 <- fgseaRes[fgseaRes$padj < 0.25,]
g2 <- g2[order(g2$NES, decreasing = T),]

g2$yax <- ifelse(g2$NES > 0, -0.02, 0.02)
g2$col <- ifelse(g2$NES > 0, "blue","red")

g2$leadingEdge <- sapply(g2$leadingEdge, toString)
g2 <- data.frame(g2)
write.table(g2, "/work/SCCC/s228508/grey_Mu_lab/choushi/gsea2.csv", sep =",", row.names = F)

# check the top2 
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

num1=6
# "#84a59d", "#f28482"
# pdf("/work/SCCC/s228508/choushi/figs/GSEA.pdf", width=5, height=4)
gseaplot2(egmt, geneSetID = "T cell activation",
          title = "adaptive immune response",
          color = "#84a59d",
          base_size = 10,
          rel_heights = c(0.8,0.2,0.2),
          subplots = 1:3,
          pvalue_table = F,
          ES_geom = "line" #line or dot
)
# dev.off()

