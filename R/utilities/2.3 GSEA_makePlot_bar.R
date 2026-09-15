
library(ggplot2)
library(dplyr)


data <- read.table("/work/SCCC/s228508/7.0 Myeloid/selected_gsea/selected_TAM_c06_inflammtory.csv", sep = ",", header = T)

data <- gsea_result2
data$yax <- ifelse(data$NES > 0, -0.02, 0.02)
data$col <- ifelse(data$NES > 0, "blue","red")


# Split data into up and downregulated pathways
data <- data %>%
  mutate(Direction = ifelse(NES > 0, "Up", "Down")) %>%
  arrange(Direction, NES)
#----
data <- g2[,c(1,6,3)]
colnames(data)<- c("pathway", "NES","FDR")
# select <- c("HALLMARK_G2M_CHECKPOINT",
#             "KEGG_CELL_CYCLE",
#             "KEGG_DNA_REPLICATION",
#             "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
#             "KEGG_P53_SIGNALING_PATHWAY",
#             "KEGG_PROSTATE_CANCER",
#             "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
#             "KEGG_HISTIDINE_METABOLISM",
#             "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
#             "HALLMARK_ANDROGEN_RESPONSE")
select <- c("HALLMARK_G2M_CHECKPOINT",
            "KEGG_CELL_CYCLE",
            "KEGG_DNA_REPLICATION",
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            # "KEGG_P53_SIGNALING_PATHWAY",
            "KEGG_PROSTATE_CANCER",
            "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
            "KEGG_HISTIDINE_METABOLISM",
            # "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
            "HALLMARK_ANDROGEN_RESPONSE")
data2 <- data[data$pathway %in% select, ]

data2 <- data2[order(data2$NES, decreasing = F),]
data2$pathway <- factor(data2$pathway, levels = data2$pathway)

figure_font_size <- 12

p <- ggplot(data2, aes(y = NES, x = reorder(pathway, NES), label = pathway)) +
  geom_bar(stat = "identity", aes(fill = NES > 0), width = 0.8) +
  geom_text(
    aes(hjust = ifelse(NES > 0, 1.05, -0.05)),  # 调整在 bar 内部
    size = figure_font_size / (14/4)
  ) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  labs(y = "Normalized Enrichment Score", x = "", title = "") +
  scale_fill_manual(
    values = c("FALSE" = "#375E97", "TRUE" = "#FB6542"), 
    labels = c("Down", "Up"), 
    name = "Expression"
  ) +
  scale_x_discrete(breaks = NULL) +
  theme_classic(base_size = figure_font_size) + 
  theme(
    panel.grid.major = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = unit(c(0.1, 0, 0.05, 0), "lines")
  )

p