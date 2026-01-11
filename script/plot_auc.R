library(tidyverse)
library(reshape2)
library(pROC)
library(reshape2)
library(cowplot)
library(vegan)

setwd("/home/dongbiao/GCN/results")
fold <- c("fold_1", "fold_2", "fold_3", "fold_4", "fold_5")
## synthetic data
metadata <- read.delim("/home/dongbiao/software/Phylo-Spec/data/Synthetic_Dataset_1/metadata.txt", check.names = FALSE, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample_id
# metadata$group <- as.factor(metadata$group) 
metadata[metadata$group == 0, "group"] <- c("Group S1 + S3")
metadata[metadata$group == 1, "group"] <- "Group S2 + S4"
n <- 1
p <- list()
for (method in c("bray_curtis", "weighted_unifrac")){
  dist <- read_qza(paste0("/home/dongbiao/software/Phylo-Spec/data/Synthetic_Dataset_1/beta/", method, "_distance_matrix.qza"))$data
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

### RF
df <- read.csv("../data/synthetic_data/results/RF_results.csv")
df[, "label"] <- metadata[df[, 1], "group"] 
RF <- c()
for (i in unique(df$fold)){
    temp <- df %>% filter(fold == i)
    roc_result <- roc(temp$label, temp$pred_prob)
    RF <- c(RF, auc(roc_result))
}
### DeepPhylo
DeepPhylo <- c()
for (i in c(1:5)){
    df <- read.csv(paste0("../data/synthetic_data/results/predictions_deepphylo_", i, ".csv"))
    df[, "label"] <- metadata[df[, 1], "group"]
    roc_result <- roc(df$label, df$Prob_1)
    DeepPhylo <- c(DeepPhylo, auc(roc_result))
}
### Phylo_spec
Phylo_Spec <- c()
for (i in c(1:5)){
    df <- read.csv(paste0("../data/synthetic_data/output/predictions_", i, ".csv"))
    df[, "label"] <- metadata[df[, 1], "group"]
    roc_result <- roc(df$label, df$Prob_Class_1)
    Phylo_Spec <- c(Phylo_Spec, auc(roc_result))
}
### GCN
GCN <- c()
for (i in c(1:5)){
    df <- read.csv(paste0("../data/synthetic_data/results/predictions_", i, ".csv"))
    roc_result <- roc(df$True_Label, df$Prob_Class_1)
    GCN <- c(GCN, auc(roc_result))
}
plot_df <- data.frame(fold = fold, RF = RF, DeepPhylo = DeepPhylo, 
                      Phylo_Spec = Phylo_Spec, PhyloGCNE = GCN)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("PhyloGCNE", "Phylo_Spec", "DeepPhylo", "RF"))

p <- ggplot(plot_df, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_line(aes(group=fold), color = "black") +
  geom_point(size = 3) + 
  labs(x = "", y = "AUC", title = "Synthetic data") +
  theme_bw() +
  theme(axis.text = element_text(size = 12))
p <- plot_grid(combined, p, align="hv", labels = c("a", "b"),
               nrow = 2, ncol=1, plot=FALSE)
ggsave(p, filename = "synthetic_results.pdf",width = 6, height = 7, useDingbats=FALSE)

## IBD 16S
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_16S_IBD//metadata.tsv", 
                     sep = "\t", row.names = 1)
### RF
df <- read.csv("../data/IBD_16S/results/RF_results.csv")
df[, "label"] <- metadata[df[, 1], "group"]
RF <- c()
for (i in unique(df$fold)){
  temp <- df %>% filter(fold == i)
  roc_result <- roc(temp$label, temp$pred_prob)
  RF <- c(RF, auc(roc_result))
}
### DeepPhylo
DeepPhylo <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/IBD_16S/results/predictions_deepphylo_", i, ".csv"))
  df[, "label"] <- metadata[df[, 1], "group"]
  roc_result <- roc(df$label, df$Prob_1)
  DeepPhylo <- c(DeepPhylo, auc(roc_result))
}
### Phylo_spec
Phylo_Spec <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/IBD_16S/output/predictions_", i, ".csv"))
  df[, "label"] <- metadata[df[, 1], "group"]
  roc_result <- roc(df$label, df$Prob_Class_1)
  Phylo_Spec <- c(Phylo_Spec, auc(roc_result))
}
### GCN
GCN <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/IBD_16S/results/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN <- c(GCN, auc(roc_result))
}
plot_df <- data.frame(fold = fold, RF = RF, DeepPhylo = DeepPhylo, 
                      Phylo_Spec = Phylo_Spec, PhyloGCNE = GCN)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("PhyloGCNE", "Phylo_Spec", "DeepPhylo", "RF"))

p1 <- ggplot(plot_df, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_line(aes(group=fold)) +
  geom_point(size = 3) + 
  labs(x = "", y = "AUC", title = "IBD 16S") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

## CRC 16S
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_16S_CRC/metadata.tsv", 
                     sep = "\t", row.names = 1)
### RF
df <- read.csv("../data/CRC_16S/results/RF_results.csv")
df[, "label"] <- metadata[df[, 1], "group"]
RF <- c()
for (i in unique(df$fold)){
  temp <- df %>% filter(fold == i)
  roc_result <- roc(temp$label, temp$pred_prob)
  RF <- c(RF, auc(roc_result))
}
### DeepPhylo
DeepPhylo <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/CRC_16S/results/predictions_deepphylo_", i, ".csv"))
  df[, "label"] <- metadata[df[, 1], "group"]
  roc_result <- roc(df$label, df$Prob_1)
  DeepPhylo <- c(DeepPhylo, auc(roc_result))
}
### Phylo_spec
Phylo_Spec <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/CRC_16S/output/predictions_", i, ".csv"))
  df[, "label"] <- metadata[df[, 1], "group"]
  roc_result <- roc(df$label, df$Prob_Class_1)
  Phylo_Spec <- c(Phylo_Spec, auc(roc_result))
}
### GCN
GCN <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/CRC_16S/results/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN <- c(GCN, auc(roc_result))
}
plot_df <- data.frame(fold = fold, RF = RF, DeepPhylo = DeepPhylo, 
                      Phylo_Spec = Phylo_Spec, PhyloGCNE = GCN)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("PhyloGCNE", "Phylo_Spec", "DeepPhylo", "RF"))

p2 <- ggplot(plot_df, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_line(aes(group=fold)) +
  geom_point(size = 3) + 
  labs(x = "", y = "AUC", title = "CRC 16S") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

### Dietary fiber 16S
metadata <- read.csv("/home/dongbiao/GCN/data/dietary_fiber/metadata.tsv", sep = "\t", row.names = 1)
### RF
df <- read.csv("../data/dietary_fiber/results/RF_results.csv")
df[, "label"] <- metadata[df[, 1], "group"]
RF <- c()
for (i in unique(df$fold)){
  temp <- df %>% filter(fold == i)
  roc_result <- roc(temp$label, temp$pred_prob)
  RF <- c(RF, auc(roc_result))
}
### DeepPhylo
DeepPhylo <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/dietary_fiber/results/predictions_deepphylo_", i, ".csv"))
  df[, "label"] <- metadata[df[, 1], "group"]
  roc_result <- roc(df$label, df$Prob_1)
  DeepPhylo <- c(DeepPhylo, auc(roc_result))
}
### Phylo_spec
Phylo_Spec <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/dietary_fiber/output/predictions_", i, ".csv"))
  df[, "label"] <- metadata[df[, 1], "group"]
  roc_result <- roc(df$label, df$Prob_Class_1)
  Phylo_Spec <- c(Phylo_Spec, auc(roc_result))
}
### GCN
GCN <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/dietary_fiber/results/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN <- c(GCN, auc(roc_result))
}
plot_df <- data.frame(fold = fold, RF = RF, DeepPhylo = DeepPhylo, 
                      Phylo_Spec = Phylo_Spec, PhyloGCNE = GCN)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("PhyloGCNE", "Phylo_Spec", "DeepPhylo", "RF"))

p3 <- ggplot(plot_df, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_line(aes(group=fold)) +
  geom_point(size = 3) + 
  labs(x = "", y = "AUC", title = "Fiber 16S") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

### CRC WGS
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_WGS_CRC/metadata.tsv", 
                     sep = "\t", row.names = 1)
### RF
df <- read.csv("../data/CRC_WGS/results/RF_results.csv")
df[, "label"] <- metadata[df[, 1], "group"]
RF <- c()
for (i in unique(df$fold)){
  temp <- df %>% filter(fold == i)
  roc_result <- roc(temp$label, temp$pred_prob)
  RF <- c(RF, auc(roc_result))
}
### DeepPhylo
DeepPhylo <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/CRC_WGS/results/predictions_deepphylo_", i, ".csv"))
  df[, "label"] <- metadata[df[, 1], "group"]
  roc_result <- roc(df$label, df$Prob_1)
  DeepPhylo <- c(DeepPhylo, auc(roc_result))
}
### Phylo_spec
Phylo_Spec <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/CRC_WGS/output/predictions_", i, ".csv"))
  df[, "label"] <- metadata[df[, 1], "group"]
  roc_result <- roc(df$label, df$Prob_Class_1)
  Phylo_Spec <- c(Phylo_Spec, auc(roc_result))
}
### GCN
GCN <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/CRC_WGS/results/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN <- c(GCN, auc(roc_result))
}
plot_df <- data.frame(fold = fold, RF = RF, DeepPhylo = DeepPhylo, 
                      Phylo_Spec = Phylo_Spec, PhyloGCNE = GCN)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("PhyloGCNE", "Phylo_Spec", "DeepPhylo", "RF"))

p4 <- ggplot(plot_df, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_line(aes(group=fold)) +
  geom_point(size = 3) + 
  labs(x = "", y = "AUC", title = "CRC WGS") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

## T2D WGS
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_WGS_T2D/metadata.tsv", 
                     sep = "\t", row.names = 1)
### RF
df <- read.csv("../data/T2D_WGS/results/RF_results.csv")
df[, "label"] <- metadata[df[, 1], "group"]
RF <- c()
for (i in unique(df$fold)){
  temp <- df %>% filter(fold == i)
  roc_result <- roc(temp$label, temp$pred_prob)
  RF <- c(RF, auc(roc_result))
}
### DeepPhylo
DeepPhylo <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/T2D_WGS/results/predictions_deepphylo_", i, ".csv"))
  df[, "label"] <- metadata[df[, 1], "group"]
  roc_result <- roc(df$label, df$Prob_1)
  DeepPhylo <- c(DeepPhylo, auc(roc_result))
}
### Phylo_spec
Phylo_Spec <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/T2D_WGS/output/predictions_", i, ".csv"))
  df[, "label"] <- metadata[df[, 1], "group"]
  roc_result <- roc(df$label, df$Prob_Class_1)
  Phylo_Spec <- c(Phylo_Spec, auc(roc_result))
}
### GCN
GCN <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/T2D_WGS/results/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN <- c(GCN, auc(roc_result))
}
plot_df <- data.frame(fold = fold, RF = RF, DeepPhylo = DeepPhylo, 
                      Phylo_Spec = Phylo_Spec, PhyloGCNE = GCN)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("PhyloGCNE", "Phylo_Spec", "DeepPhylo", "RF"))

p5 <- ggplot(plot_df, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_line(aes(group=fold)) +
  geom_point(size = 3) + 
  labs(x = "", y = "AUC", title = "T2D WGS") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

## Multi-classification
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_Multi-classification/metadata.tsv", 
                     sep = "\t", row.names = 1)
metadata_unique <- metadata[!duplicated(metadata$disease), ]
### RF
df <- read.csv("../data/Multi_classification/results/RF_results.csv")
colnames(df) <- c("sample_id", "fold", "ASD", "CRC", "HC", "IBD", "IBS")
df <- df %>% melt(id.vars = c("sample_id", "fold"))
RF <- c()
for (j in c("ASD", "CRC", "HC", "IBD", "IBS")){
  for (i in unique(df$fold)){
    temp <- df %>% filter(fold == i) %>% filter(variable == j)
    temp[, "label"] <- 0
    temp[metadata[temp[, 1], "disease"] == j, "label"] <- 1
    roc_result <- roc(temp$label, temp$value)
    RF <- c(RF, auc(roc_result))
  }
}

### DeepPhylo
DeepPhylo <- c()
for (j in c("ASD", "CRC", "HC", "IBD", "IBS")){
for (i in c(1:5)){
  df <- read.csv(paste0("../data/Multi_classification/results/predictions_deepphylo_", i, ".csv"))
  colnames(df) <- c("sample_id", "ASD", "CRC", "HC", "IBD", "IBS")
  df <- df %>% melt(id.vars = c("sample_id")) %>% filter(variable == j)
  df[, "label"] <- 0
  df[metadata[df[, "sample_id"], "disease"] == j, "label"] <- 1
  roc_result <- roc(df$label, df$value)
  DeepPhylo <- c(DeepPhylo, auc(roc_result))
}
}
### Phylo_spec
Phylo_Spec <- c()
for (j in c("ASD", "CRC", "HC", "IBD", "IBS")){
for (i in c(1:5)){
  df <- read.csv(paste0("../data/Multi_classification/output/predictions_", i, ".csv"))
  colnames(df) <- c("sample_id", "ASD", "CRC", "HC", "IBD", "IBS")
  df <- df %>% melt(id.vars = c("sample_id")) %>% filter(variable == j)
  df[, "label"] <- 0
  df[metadata[df[, "sample_id"], "disease"] == j, "label"] <- 1
  roc_result <- roc(df$label, df$value)
  Phylo_Spec <- c(Phylo_Spec, auc(roc_result))
}
}
### GCN
GCN <- c()
for (j in c("ASD", "CRC", "HC", "IBD", "IBS")){
for (i in c(1:5)){
  df <- read.csv(paste0("../data/Multi_classification/results/predictions_", i, ".csv"))
  colnames(df) <- c("sample_id", "ASD", "CRC", "HC", "IBD", "IBS", "label")
  df <- df %>% melt(id.vars = c("sample_id", "label")) %>% filter(variable == j)
  df[, "label"] <- 0
  df[metadata[df[, "sample_id"], "disease"] == j, "label"] <- 1
  roc_result <- roc(df$label, df$value)
  GCN <- c(GCN, auc(roc_result))
}}

plot_df <- data.frame(fold = fold, RF = RF, DeepPhylo = DeepPhylo, 
                      Phylo_Spec = Phylo_Spec, PhyloGCNE = GCN)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("PhyloGCNE", "Phylo_Spec", "DeepPhylo", "RF"))
plot_df <- plot_df %>% group_by(fold, variable) %>% summarise(mean_value = mean(value))
p6 <- ggplot(plot_df, aes(x = variable, y = mean_value)) +
  geom_boxplot() +
  geom_line(aes(group=fold))+
  geom_point(size = 3) + 
  labs(x = "", y = "AUC", title = "Multiclass") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

p <- plot_grid(p1, p2, p3, p4, p5, p6, align="hv", labels = c("a", "b", "c", "d", "e", "f"),
               nrow = 2, ncol=3, plot=FALSE)
ggsave(p, filename = "Benchmark_auc.pdf",width = 10, height = 7, useDingbats=FALSE)

