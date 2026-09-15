rm(list = ls())

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(Seurat)
library(DESeq2)
library(dplyr)
library(tibble)
library(ggrepel)
library(apeglm)
library(limma)

# seurat.harmony <- readRDS("/work/SCCC/s228508/grey_Mu_lab/choushi/primary_harmony(240929).rds")
# seurat.harmony <- readRDS("/work/SCCC/s228508/grey_Mu_lab/choushi/epi_rough(240930).rds")
tumor.epi <- readRDS("/gpfs/gibbs/project/mu/qx75/Grey/grey_Mu_lab/choushi/epi_rough(241023).rds")
DefaultAssay(tumor.epi) <- "RNA" 
# Idents(seurat.harmony) <- "Sample_ID"
Idents(tumor.epi) <- "epi_rough"

# markers <- FindAllMarkers(object = seurat.harmony,
#                              only.pos = T,
#                              min.pct = 0.1,
#                              logfc.threshold = 0.01)
# markers = subset(ann_markers, ann_markers$cluster %in% c("2"))


markers <-FindMarkers(object = tumor.epi,
                      ident.1 = "c2", #ident.2 = "sgNT",
                      only.pos = F,
                      min.pct = 0.1,
                      logfc.threshold = 0.5)
up_deg <- rownames(markers)

#--------------
markers_2 <- markers %>%
  mutate(Type = if_else(p_val_adj > 0.05,
                        "ns",
                        if_else(abs(avg_log2FC) <= 0.5, "ns",
                                if_else(avg_log2FC > 0.5, "up", "down")))) %>%
  arrange(desc(abs(avg_log2FC))) %>% rownames_to_column("Gene_Symbol")

up_deg <- subset(markers_2, markers_2$Type %in% "up")
length(rownames(up_deg))



#----------------
degs.list=up_deg$Gene_Symbol
degs.list
erich.go.BP = enrichGO(gene = degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "ALL",
                       pAdjustMethod = 'fdr',
                       pvalueCutoff = 1,
                       qvalueCutoff = 1)

result_go=erich.go.BP@result 
result_go2 <- subset(result_go, p.adjust < 0.05)
result_go2 <- result_go2[order(result_go2$p.adjust, decreasing = FALSE),]

dotplot(erich.go.BP,showCategory = 20)
barplot(erich.go.BP,showCategory = 30)

# Order by p-value (ascending)
top20_result_go <- result_go[order(result_go$pvalue), ][1:35, ]

# Create a new enrichGO object with only the top 20 results
top20_erich.go.BP <- erich.go.BP
top20_erich.go.BP@result <- top20_result_go

# Plot using dotplot and barplot
dotplot(top20_erich.go.BP, showCategory = 20)
barplot(top20_erich.go.BP, showCategory = 20)

#-----------
# Extract and order the results by p-value
result_go <- erich.go.BP@result
top20_result_go <- result_go2 %>%
  arrange(qvalue) %>%
  slice(1:40)

# Prepare the data frame for ggplot
top20_result_go <- top20_result_go %>%
  mutate(Description = factor(Description, levels = rev(Description)))  # Reverse order for better plotting

# Plot with ggplot2 (barplot example)
ggplot(top20_result_go, aes(x = Description, y = -log10(qvalue))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip the axes for better readability
  labs(x = "GO Biological Process", y = "-log10(q-value)", title = "Top 40 GO Biological Processes") +
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 12),       # Adjust x-axis text size
    axis.text.y = element_text(size = 12),       # Adjust y-axis text size
    axis.title.x = element_text(size = 14),      # Adjust x-axis title size
    axis.title.y = element_text(size = 14),      # Adjust y-axis title size
    plot.title = element_text(size = 16, face = "bold")  # Adjust title size
  )

# Plot with ggplot2 (dotplot example)
ggplot(top20_result_go, aes(x = -log10(qvalue), y = Description)) +
  geom_point(aes(size = Count, color = -log10(qvalue))) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "-log10(q-value)", y = "GO Biological Process", title = "Top 20 GO Biological Processes") +
  theme_minimal()





#-------
erich.go.MF = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
result_mf=erich.go.MF@result 
dotplot(erich.go.MF,showCategory = 20)
barplot(erich.go.MF,showCategory = 20)

# Extract and order the results by p-value
result_mf <- erich.go.MF@result
top20_result_mf <- result_mf %>%
  arrange(qvalue) %>%
  slice(1:20)

# Prepare the data frame for ggplot
top20_result_mf <- top20_result_mf %>%
  mutate(Description = factor(Description, levels = rev(Description)))  # Reverse order for better plotting

# Plot with ggplot2 (barplot example)
ggplot(top20_result_mf, aes(x = Description, y = -log10(qvalue))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip the axes for better readability
  labs(x = "GO MF", y = "-log10(q-value)", title = "Top 20 MF") +
  theme_minimal()

#---------
erich.go.CC = enrichGO(gene =degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
dotplot(erich.go.CC,showCategory = 20)
barplot(erich.go.CC,showCategory = 20)

# Extract and order the results by p-value
result_cc <- erich.go.CC@result
top20_result_cc <- result_cc %>%
  arrange(pvalue) %>%
  slice(1:20)

# Prepare the data frame for ggplot
top20_result_cc <- top20_result_cc %>%
  mutate(Description = factor(Description, levels = rev(Description)))  # Reverse order for better plotting

# Plot with ggplot2 (barplot example)
ggplot(top20_result_cc, aes(x = Description, y = -log10(qvalue))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip the axes for better readability
  labs(x = "GO CC", y = "-log10(q-value)", title = "Top 20 CC") +
  theme_minimal()


#---------
#KEGG
#ID change：
degs.list=up_deg$Gene_Symbol
degs.list=up_deg$gene
degs.list
# degs.list <- up_deg
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys =  degs.list,
                       keytype = "SYMBOL",
                       column = "ENTREZID")

enrich.kegg.res <- enrichKEGG(gene = DEG.entrez_id,
                              organism = "hsa", # mmu
                              keyType = "kegg")
dotplot(enrich.kegg.res)

# library(R.utils)
# R.utils::setOption('clusterProfiler.download.method','auto')

result_kegg <- enrich.kegg.res@result

write.table(result_kegg, "enrich_kegg.csv",sep = ",",col.names = T, row.names = T)





