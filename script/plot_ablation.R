library(tidyverse)
library(reshape2)
library(pROC)
library(reshape2)
library(cowplot)

setwd("/home/dongbiao/GCN/results")
fold <- c("fold_1", "fold_2", "fold_3", "fold_4", "fold_5")
## IBD 16S
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_16S_IBD//metadata.tsv", 
                     sep = "\t", row.names = 1)
### baseline
### GCN
GCN <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/IBD_16S/results/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  Phylo_Spec <- c(Phylo_Spec, auc(roc_result))
  GCN <- c(GCN, auc(roc_result))
}
### GCN ablation distance
GCN_distance <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/IBD_16S/ablation_distance/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN_distance <- c(GCN_distance, auc(roc_result))
}
### GCN ablation phylo
GCN_phylo <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/IBD_16S/ablation_phylo/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN_phylo <- c(GCN_phylo, auc(roc_result))
}

plot_df <- data.frame(fold = fold, Full_tree = GCN, Topology_Only = GCN_distance, 
                      Tips_only = GCN_phylo)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("Full_tree", "Topology_Only", "Tips_only"))
plot_df <- plot_df %>% filter(variable != "Topology_Only")

p1 <- ggplot(plot_df, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_line(aes(group=fold)) +
  geom_point(size = 3) + 
  labs(x = "", y = "AUC", title = "IBD_16S") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

## CRC 16S
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_16S_CRC/metadata.tsv", 
                     sep = "\t", row.names = 1)
### baseline
### GCN
GCN <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/CRC_16S/results/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  Phylo_Spec <- c(Phylo_Spec, auc(roc_result))
  GCN <- c(GCN, auc(roc_result))
}
### GCN ablation distance
GCN_distance <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/CRC_16S/ablation_distance/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN_distance <- c(GCN_distance, auc(roc_result))
}
### GCN ablation phylo
GCN_phylo <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/CRC_16S/ablation_phylo/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN_phylo <- c(GCN_phylo, auc(roc_result))
}

plot_df <- data.frame(fold = fold, Full_tree = GCN, Topology_Only = GCN_distance, 
                      Tips_only = GCN_phylo)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("Full_tree", "Topology_Only", "Tips_only"))
plot_df <- plot_df %>% filter(variable != "Topology_Only")
p2 <- ggplot(plot_df, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_line(aes(group=fold)) +
  geom_point(size = 3) + 
  labs(x = "", y = "AUC", title = "CRC_16S") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))


## fiber 16S
metadata <- read.csv("/home/dongbiao/GCN/data/dietary_fiber/metadata.tsv", 
                     sep = "\t", row.names = 1)
### baseline
### GCN
GCN <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/dietary_fiber/results/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  Phylo_Spec <- c(Phylo_Spec, auc(roc_result))
  GCN <- c(GCN, auc(roc_result))
}
### GCN ablation distance
GCN_distance <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/dietary_fiber/ablation_distance/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN_distance <- c(GCN_distance, auc(roc_result))
}
### GCN ablation phylo
GCN_phylo <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/dietary_fiber/ablation_phylo/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN_phylo <- c(GCN_phylo, auc(roc_result))
}

plot_df <- data.frame(fold = fold, Full_tree = GCN, Topology_Only = GCN_distance, 
                      Tips_only = GCN_phylo)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("Full_tree", "Topology_Only", "Tips_only"))
plot_df <- plot_df %>% filter(variable != "Topology_Only")

p3 <- ggplot(plot_df, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_line(aes(group=fold)) +
  geom_point(size = 3) + 
  labs(x = "", y = "AUC", title = "Fiber_16S") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

### CRC WGS
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_WGS_CRC/metadata.tsv", 
                     sep = "\t", row.names = 1)
### baseline
### GCN
GCN <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/CRC_WGS/results/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  Phylo_Spec <- c(Phylo_Spec, auc(roc_result))
  GCN <- c(GCN, auc(roc_result))
}
### GCN ablation distance
GCN_distance <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/CRC_WGS/ablation_distance/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN_distance <- c(GCN_distance, auc(roc_result))
}
### GCN ablation phylo
GCN_phylo <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/CRC_WGS/ablation_phylo/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN_phylo <- c(GCN_phylo, auc(roc_result))
}

plot_df <- data.frame(fold = fold, Full_tree = GCN, Topology_Only = GCN_distance, 
                      Tips_only = GCN_phylo)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("Full_tree", "Topology_Only", "Tips_only"))
plot_df <- plot_df %>% filter(variable != "Topology_Only")
p4 <- ggplot(plot_df, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_line(aes(group=fold)) +
  geom_point(size = 3) + 
  labs(x = "", y = "AUC", title = "CRC_WGS") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

### T2D WGS
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_WGS_T2D/metadata.tsv", 
                     sep = "\t", row.names = 1)
### baseline
### GCN
GCN <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/T2D_WGS/results/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  Phylo_Spec <- c(Phylo_Spec, auc(roc_result))
  GCN <- c(GCN, auc(roc_result))
}
### GCN ablation distance
GCN_distance <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/T2D_WGS/ablation_distance/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN_distance <- c(GCN_distance, auc(roc_result))
}
### GCN ablation phylo
GCN_phylo <- c()
for (i in c(1:5)){
  df <- read.csv(paste0("../data/T2D_WGS/ablation_phylo/predictions_", i, ".csv"))
  roc_result <- roc(df$True_Label, df$Prob_Class_1)
  GCN_phylo <- c(GCN_phylo, auc(roc_result))
}

plot_df <- data.frame(fold = fold, Full_tree = GCN, Topology_Only = GCN_distance, 
                      Tips_only = GCN_phylo)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("Full_tree", "Topology_Only", "Tips_only"))
plot_df <- plot_df %>% filter(variable != "Topology_Only")
p5 <- ggplot(plot_df, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_line(aes(group=fold)) +
  geom_point(size = 3) + 
  labs(x = "", y = "AUC", title = "CRC_WGS") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))


### Multi-classification
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_Multi-classification/metadata.tsv", 
                     sep = "\t", row.names = 1)
### baseline
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
### GCN ablation distance
GCN_distance <- c()
for (j in c("ASD", "CRC", "HC", "IBD", "IBS")){
  for (i in c(1:5)){
    df <- read.csv(paste0("../data/Multi_classification/ablation_distance/predictions_", i, ".csv"))
    colnames(df) <- c("sample_id", "ASD", "CRC", "HC", "IBD", "IBS", "label")
    df <- df %>% melt(id.vars = c("sample_id", "label")) %>% filter(variable == j)
    df[, "label"] <- 0
    df[metadata[df[, "sample_id"], "disease"] == j, "label"] <- 1
    roc_result <- roc(df$label, df$value)
    GCN_distance <- c(GCN_distance, auc(roc_result))
  }}
### GCN ablation phylo
GCN_phylo <- c()
for (j in c("ASD", "CRC", "HC", "IBD", "IBS")){
  for (i in c(1:5)){
    df <- read.csv(paste0("../data/Multi_classification/ablation_phylo/predictions_", i, ".csv"))
    colnames(df) <- c("sample_id", "ASD", "CRC", "HC", "IBD", "IBS", "label")
    df <- df %>% melt(id.vars = c("sample_id", "label")) %>% filter(variable == j)
    df[, "label"] <- 0
    df[metadata[df[, "sample_id"], "disease"] == j, "label"] <- 1
    roc_result <- roc(df$label, df$value)
    GCN_phylo <- c(GCN_phylo, auc(roc_result))
  }}

plot_df <- data.frame(fold = fold, Full_tree = GCN, Topology_Only = GCN_distance, 
                      Tips_only = GCN_phylo)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("Full_tree", "Topology_Only", "Tips_only"))
plot_df <- plot_df %>% group_by(fold, variable) %>% summarise(mean_value = mean(value))
plot_df <- plot_df %>% filter(variable != "Topology_Only")
p6 <- ggplot(plot_df, aes(x = variable, y = mean_value)) +
  geom_boxplot() +
  geom_line(aes(group=fold)) +
  geom_point(size = 3) + 
  labs(x = "", y = "AUC", title = "Multiclass") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
p <- plot_grid(p1, p2, p3, p4, p5, p6, align="hv", labels = c("a", "c", "b", "d", "e", "f"),
               nrow = 2, ncol=3, plot=FALSE)
ggsave(p, filename = "Ablation_auc.pdf",width = 10, height = 7, useDingbats=FALSE)

