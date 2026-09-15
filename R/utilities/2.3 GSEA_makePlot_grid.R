library(ggplot2)
library(dplyr)

gsea_combined <- read.table("/work/SCCC/s228508/7.0 Myeloid/combined2.csv", sep = ",", header = T)

# Assuming your dataframe is named df and column A is named "A"
gsea_combined$Significance <- ifelse(gsea_combined$p.adjust == "ns", "ns", "Significant")

gsea_combined$NES_Adjusted <- ifelse(gsea_combined$Significance == "ns", NA, gsea_combined$NES)
gsea_combined$NES_Adjusted <- as.numeric(gsea_combined$NES_Adjusted)


p <-
  ggplot(gsea_combined, aes(x = cluster, y = pathway)) +
  geom_tile(aes(fill = NES_Adjusted)) +  # Use the adjusted NES column
  scale_fill_gradient2(low = "#3B53A4", high = "#8F2125", mid = "white",  #low = "#5B76A2", high = "#8F2125"
                       midpoint = 0, limits = c(-2.5, 2.5), 
                       name = "NES (GSEA, adj-p<0.05)", 
                       na.value = "white") +  # Set white for non-significant (NA) values
  # geom_text(aes(label = ifelse(Significance == "Significant", "*", "")), 
  #           color = "black", size = 4, vjust = 1.5) +  # Asterisks for significant values
  labs(x = NULL, y = NULL) +  # Remove axis labels
  theme_minimal() +  # Use a minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank())  # Remove gridlines
p
ggsave("/work/SCCC/s228508/7.0 Myeloid/GSEA_combined/GSEA_combined.pdf", plot = p, height = 22, width = 9)


#--------
# have lines
gsea_combined2 <- subset(gsea_combined, gsea_combined$Function %in% c("others"))

ggplot(gsea_combined2, aes(x = cluster, y = pathway)) +
  geom_tile(aes(fill = NES_Adjusted), color = "grey") +  # Add black lines around each square
  scale_fill_gradient2(low = "#3B53A4", high = "#ED2C24", mid = "white", 
                       midpoint = 0, limits = c(-2.5, 2.5), 
                       name = "NES (GSEA, adj-p<0.05)", 
                       na.value = "white") +  # Set white for non-significant (NA) values
  # geom_text(aes(label = ifelse(Significance == "Significant", "*", "")), 
  #           color = "black", size = 4, vjust = 1.5) +  # Asterisks for significant values
  labs(x = NULL, y = NULL) +  # Remove axis labels
  theme_minimal() +  # Use a minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank())  # Remove gridlines

