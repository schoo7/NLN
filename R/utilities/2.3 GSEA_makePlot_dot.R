library(ggplot2)

data <- gsea_result2[,c(2,5,7)]
colnames(data)<- c("pathway", "NES","FDR")

data <- g2[,c(1,6,3)]
colnames(data)<- c("pathway", "NES","FDR")
select <- c("HALLMARK_G2M_CHECKPOINT",
            "KEGG_CELL_CYCLE",
            "KEGG_DNA_REPLICATION",
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            "KEGG_P53_SIGNALING_PATHWAY",
            "KEGG_PROSTATE_CANCER",
            "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
            "KEGG_HISTIDINE_METABOLISM",
            "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY",
            "HALLMARK_ANDROGEN_RESPONSE")
data2 <- data[data$pathway %in% select, ]

data2 <- data2[order(data2$NES, decreasing = F),]
data2$pathway <- factor(data2$pathway, levels = data2$pathway)

ggplot(data2, aes(x = NES, y = pathway)) +
  # Add plaid background with horizontal and vertical lines
  geom_hline(yintercept = seq(1, length(data2$pathway), by = 1), color = "grey90") + 
  geom_vline(xintercept = seq(-3, 3, by = 1), color = "grey90") +
  geom_point(aes(size = -log10(FDR), fill = NES > 0), shape = 21, color = "black") +
  scale_size_continuous(range = c(2, 10), breaks = c(2, 5, 7.5), name = "-Log10(FDR)") +
  scale_fill_manual(values = c("#6593CD", "#EB5D6B")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "italic"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Add a black frame
  ) +
  labs(x = "NES", y = NULL, title = "sgNLN_WNT_vs_sgNLN") +
  scale_x_continuous(limits = c(-3, 3),breaks = seq(-3, 3, by = 1))  # Specify x-axis ticks

