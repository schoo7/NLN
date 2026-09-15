library(ggplot2)

# Violin plot
ggplot(slingshot_df, aes(x = epi_rough, y = slingPseudotime_2, fill = epi_rough)) +
  geom_violin(trim = FALSE, width = 1, adjust = 1.2, alpha = 0.8, color = "black", size = 0.3) +
  geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.4) +
  theme_classic() +
  labs(x = "Clusters", y = "Pseudotime") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none")

ggplot(slingshot_df, aes(x = epi_rough, y = slingPseudotime_2, fill = epi_rough)) +
  geom_violin(trim = FALSE, width = 1, adjust = 1.2, alpha = 0.8, color = "black", size = 0.3) +
  geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.4) +
  theme_classic() +
  labs(x = "Clusters", y = "Pseudotime") +
  scale_fill_brewer(palette = "Set2") +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black")
  )

#-------
library(ggplot2)
library(RColorBrewer)
my_color <- c("#d53e4f", "#bf812d", "#5ab4ac", "#4575b4")
ggplot(slingshot_df, aes(x = epi_rough, y = slingPseudotime_2, fill = epi_rough)) +
  geom_violin(
    trim = FALSE,
    width = 1,
    adjust = 1.2,
    alpha = 0.8,
    color = "black",
    size = 0.3
  ) +
  geom_boxplot(
    width = 0.2,
    outlier.size = 0.5,
    alpha = 0.4,
    color = "black"
  ) +
  scale_fill_manual(values = my_color) +
  labs(
    x = "Clusters",
    y = "Pseudotime",
    title = "Pseudotime Distribution Across Clusters"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line = element_line(color = "black")
  )
