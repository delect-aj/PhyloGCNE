library(tidyverse)
library(reshape2)
library(pROC)
library(reshape2)
library(cowplot)
library(patchwork)

setwd("/home/dongbiao/GCN/results")
fold <- c("fold_1", "fold_2", "fold_3", "fold_4", "fold_5")
percentile <- c(0.1, 0.2, 0.4, 0.5, 0.6, 0.8)
## IBD 16S
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_16S_IBD/metadata.tsv", 
                     sep = "\t", row.names = 1)

df <- read.csv("../data/IBD_16S/results/RF_results.csv")
df[, "label"] <- metadata[df[, 1], "group"] 
RF <- c()
for (i in unique(df$fold)){
  temp <- df %>% filter(fold == i)
  roc_result <- roc(temp$label, temp$pred_prob)
  RF <- c(RF, auc(roc_result))
}
mean_RF <- mean(RF)
plot_df <- read.csv("/home/dongbiao/GCN/data/IBD_16S/results/RF_results_biomark.csv")
plot_df <- plot_df %>% mutate(group_1 = case_when(
      group_1 == "without_low" ~ "bottom-k exclusion",
      group_1 == "without_high" ~ "top-k exclusion",
      TRUE ~ group_1 
    ))
p1 <- ggplot(plot_df, aes(x = group_2, y = value, color = group_1)) +
  geom_line() +
  geom_hline(yintercept = mean_RF, linetype = "dashed", color = "black", size = 0.5) +
  geom_point(size = 3) + 
  labs(x = "Percentile", y = "AUC", title = "IBD_16S", color = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 12))

## CRC 16S
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_16S_CRC/metadata.tsv", 
                     sep = "\t", row.names = 1)
df <- read.csv("../data/CRC_16S/results/RF_results.csv")
df[, "label"] <- metadata[df[, 1], "group"] 
RF <- c()
for (i in unique(df$fold)){
  temp <- df %>% filter(fold == i)
  roc_result <- roc(temp$label, temp$pred_prob)
  RF <- c(RF, auc(roc_result))
}
mean_RF <- mean(RF)
plot_df <- read.csv("/home/dongbiao/GCN/data/CRC_16S/results/RF_results_biomark.csv")
plot_df <- plot_df %>% mutate(group_1 = case_when(
  group_1 == "without_low" ~ "bottom-k exclusion",
  group_1 == "without_high" ~ "top-k exclusion",
  TRUE ~ group_1 
))
p2 <- ggplot(plot_df, aes(x = group_2, y = value, color = group_1)) +
  geom_line() +
  geom_hline(yintercept = mean_RF, linetype = "dashed", color = "black", size = 0.5) +
  geom_point(size = 3) + 
  labs(x = "Percentile", y = "AUC", title = "CRC_16S", color = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 12))

### Fiber 16S
metadata <- read.csv("/home/dongbiao/GCN/data/dietary_fiber/metadata.tsv", 
                     sep = "\t", row.names = 1)
df <- read.csv("../data/dietary_fiber/results/RF_results.csv")
df[, "label"] <- metadata[df[, 1], "group"] 
RF <- c()
for (i in unique(df$fold)){
  temp <- df %>% filter(fold == i)
  roc_result <- roc(temp$label, temp$pred_prob)
  RF <- c(RF, auc(roc_result))
}
mean_RF <- mean(RF)
plot_df <- read.csv("/home/dongbiao/GCN/data/dietary_fiber/results/RF_results_biomark.csv")
plot_df <- plot_df %>% mutate(group_1 = case_when(
  group_1 == "without_low" ~ "bottom-k exclusion",
  group_1 == "without_high" ~ "top-k exclusion",
  TRUE ~ group_1 
))
p3 <- ggplot(plot_df, aes(x = group_2, y = value, color = group_1)) +
  geom_line() +
  geom_hline(yintercept = mean_RF, linetype = "dashed", color = "black", size = 0.5) +
  geom_point(size = 3) + 
  labs(x = "Percentile", y = "AUC", title = "Fiber_16S", color = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 12))


### CRC WGS
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_WGS_CRC/metadata.tsv", 
                     sep = "\t", row.names = 1)
df <- read.csv("../data/CRC_WGS/results/RF_results.csv")
df[, "label"] <- metadata[df[, 1], "group"] 
RF <- c()
for (i in unique(df$fold)){
  temp <- df %>% filter(fold == i)
  roc_result <- roc(temp$label, temp$pred_prob)
  RF <- c(RF, auc(roc_result))
}
mean_RF <- mean(RF)
plot_df <- read.csv("/home/dongbiao/GCN/data/CRC_WGS/results/RF_results_biomark.csv")
plot_df <- plot_df %>% mutate(group_1 = case_when(
  group_1 == "without_low" ~ "bottom-k exclusion",
  group_1 == "without_high" ~ "top-k exclusion",
  TRUE ~ group_1 
))
p4 <- ggplot(plot_df, aes(x = group_2, y = value, color = group_1)) +
  geom_line() +
  geom_hline(yintercept = mean_RF, linetype = "dashed", color = "black", size = 0.5) +
  geom_point(size = 3) + 
  labs(x = "Percentile", y = "AUC", title = "CRC_WGS", color = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 12))

### T2D WGS
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_WGS_T2D/metadata.tsv", 
                     sep = "\t", row.names = 1)
df <- read.csv("../data/T2D_WGS/results/RF_results.csv")
df[, "label"] <- metadata[df[, 1], "group"] 
RF <- c()
for (i in unique(df$fold)){
  temp <- df %>% filter(fold == i)
  roc_result <- roc(temp$label, temp$pred_prob)
  RF <- c(RF, auc(roc_result))
}
mean_RF <- mean(RF)
plot_df <- read.csv("/home/dongbiao/GCN/data/T2D_WGS/results/RF_results_biomark.csv")
plot_df <- plot_df %>% mutate(group_1 = case_when(
  group_1 == "without_low" ~ "bottom-k exclusion",
  group_1 == "without_high" ~ "top-k exclusion",
  TRUE ~ group_1 
))
p5 <- ggplot(plot_df, aes(x = group_2, y = value, color = group_1)) +
  geom_line() +
  geom_hline(yintercept = mean_RF, linetype = "dashed", color = "black", size = 0.5) +
  geom_point(size = 3) + 
  labs(x = "Percentile", y = "AUC", title = "T2D_WGS", color = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 12))

### Mult-Classification
metadata <- read.csv("/home/dongbiao/software/Phylo-Spec/data/Real_Dateset_Multi-classification/metadata.tsv", 
                     sep = "\t", row.names = 1)

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
mean_RF <- mean(RF)
plot_df <- read.csv("/home/dongbiao/GCN/data/Multi_classification/results/RF_results_biomark.csv")
plot_df <- plot_df %>% mutate(group_1 = case_when(
  group_1 == "without_low" ~ "bottom-k exclusion",
  group_1 == "without_high" ~ "top-k exclusion",
  TRUE ~ group_1 
))
p6 <- ggplot(plot_df, aes(x = group_2, y = value, color = group_1)) +
  geom_line() +
  geom_hline(yintercept = mean_RF, linetype = "dashed", color = "black", size = 0.5) +
  geom_point(size = 3) + 
  labs(x = "Percentile", y = "AUC", title = "Multiclass", color = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 12))

p <- (p1 + p2 + p3 + p4 + p5 + p6) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom') 
p <- p + plot_annotation(tag_levels = list(c('a', 'b', 'c', 'd', 'e', 'f')))
ggsave(p, filename = "Biomarker_auc.pdf", width = 10, height = 7, useDingbats=FALSE)

