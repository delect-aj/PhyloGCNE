library(tidyverse)
library(reshape2)


fold <- c("fold_1", "fold_2", "fold_3", "fold_4", "fold_5")
RF <- c(0.84, 0.81, 0.78, 0.79, 0.73)
DeepPhylo <- c(0.84, 0.86, 0.66, 0.94, 0.94)
Phylo_Spec <- c(0.94, 0.94, 0.93, 0.91, 0.94)
GCN <- c(1, 0.97, 0.95, 1, 0.99)

plot_df <- data.frame(fold = fold, RF = RF, DeepPhylo = DeepPhylo, 
                      Phylo_Spec = Phylo_Spec, GCN = GCN)
plot_df <- plot_df %>% melt(id.vars = "fold")
plot_df$variable <- factor(plot_df$variable, levels = c("GCN", "Phylo_Spec", "DeepPhylo", "RF"))

ggplot(plot_df, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_point(size = 3) + 
  labs(x = "", y = "AUC") +
  theme_bw() +
  theme(axis.text = element_text(size = 12))
