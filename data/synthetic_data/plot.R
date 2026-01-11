library(tidyverse)
library(vegan)
library(qiime2R)
library(patchwork)

# setwd("/home/dongbiao/GCN/data/synthetic_data/")
setwd("/home/dongbiao/software/Phylo-Spec/data/Synthetic_Dataset_1")
metadata <- read.delim("metadata.txt", check.names = FALSE, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample_id
# metadata$group <- as.factor(metadata$group) 
metadata[metadata$group == 0, "group"] <- c("Group S1 + S3")
metadata[metadata$group == 1, "group"] <- "Group S2 + S4"
n <- 1
p <- list()
for (method in c("bray_curtis", "weighted_unifrac")){
    dist <- read_qza(paste0("beta/", method, "_distance_matrix.qza"))$data
    dist <- as.dist(dist)
    adonis_result <- adonis2(dist ~ group, data = metadata, permutations = 999)
    R2_val <- round(adonis_result$R2[1], 3)
    P_val  <- adonis_result$`Pr(>F)`[1]
    stats_text <- paste("PERMANOVA:\nR2 =", R2_val, "\nP =", P_val)
    
    pcoa_bc <- cmdscale(dist, k = 2, eig = TRUE)
    var_exp <- round(pcoa_bc$eig / sum(pcoa_bc$eig) * 100, 1)
    PC1_label <- paste0("PC1 (", var_exp[1], "%)")
    PC2_label <- paste0("PC2 (", var_exp[2], "%)")
    plot_data <- data.frame(PC1 = pcoa_bc$points[,1],
                            PC2 = pcoa_bc$points[,2])
    plot_data <- merge(plot_data, metadata, by = "row.names")
    
    stats_text <- paste("R2 =", R2_val, "\n P =", P_val)
    
    p[[n]] <- ggplot(plot_data, aes(x = PC1, y = PC2, color = group, fill = group)) +
      geom_point(size = 3, alpha = 0.5) + 
      stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.97, show.legend = FALSE) +
      labs(x = PC1_label, y = PC2_label, title = method, color = NULL, fill = NULL) +
      annotate("text", 
               x = max(plot_data$PC1), 
               y = max(plot_data$PC2),
               label = stats_text, hjust = 1, vjust = 1, size = 4) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 12),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12)) +
      guides(color = guide_legend(override.aes = list(size = 3)))
    
    n <- n + 1
}

combined <- p[[1]] + p[[2]] + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("beta/beta_res.pdf", combined, width = 6, height = 3.5)
