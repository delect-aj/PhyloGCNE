library(tidyverse)
library(reshape2)
library(reshape2)
library(cowplot)
library(patchwork)
library(stringr)

setwd("/home/dongbiao/GCN/results")
### greengenes1 tax
tax_ref <- read.csv("/beegfs/db/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt", 
                    sep="\t", header = FALSE)
tax_ref$V1 <- as.character(tax_ref$V1)
taxonomy_df <- as.data.frame(str_split_fixed(tax_ref$V2, pattern = "; ", n = 7))
colnames(taxonomy_df) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(taxonomy_df) <- tax_ref$V1

### greengenes2 tax
tax_ref <- read.csv("/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv", 
                    sep="\t")
taxonomy_df_2 <- as.data.frame(str_split_fixed(tax_ref$Taxon, pattern = "; ", n = 7))
colnames(taxonomy_df_2) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(taxonomy_df_2) <- tax_ref$Feature.ID

## IBD 16S
importance_fid <- read.csv("../data/IBD_16S/results/importance_fid.csv") %>% 
  arrange(importance_score) %>% tail(10)
importance_fid$X <- as.character(importance_fid$X)
colnames(importance_fid)[1] <- "fid"
importance_fid <- cbind(importance_fid, taxonomy_df[importance_fid$fid, ])
importance_fid$fid <- factor(importance_fid$fid, levels = importance_fid$fid)
importance_fid <- importance_fid %>%
  mutate(
    tax = case_when(
      Species != "s__" ~ paste0(Species),
      Genus != "g__" ~ paste0(Genus),
      Family != "f__" ~ paste0(Family),
      Order != "o__" ~ paste0(Order),
      Class != "c__" ~ paste0(Class),
      Phylum != "p__" ~ paste0(Phylum),
      TRUE ~ "unclassified")
  )
labels_mapping <- importance_fid$tax
names(labels_mapping) <- importance_fid$fid
p1 <- ggplot(importance_fid, aes(x = importance_score, y = fid)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "IBD_16S", x = "Importance_score", y = "") +
  scale_y_discrete(labels = labels_mapping) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12))


### CRC 16S
importance_fid <- read.csv("../data/CRC_16S/results/importance_fid.csv") %>% 
  arrange(importance_score) %>% tail(10)
importance_fid$X <- as.character(importance_fid$X)
colnames(importance_fid)[1] <- "fid"
importance_fid <- cbind(importance_fid, taxonomy_df[importance_fid$fid, ])
importance_fid$fid <- factor(importance_fid$fid, levels = importance_fid$fid)
importance_fid <- importance_fid %>%
  mutate(
    tax = case_when(
      Species != "s__" ~ paste0(Species),
      Genus != "g__" ~ paste0(Genus),
      Family != "f__" ~ paste0(Family),
      Order != "o__" ~ paste0(Order),
      Class != "c__" ~ paste0(Class),
      Phylum != "p__" ~ paste0(Phylum),
      TRUE ~ "unclassified")
  )
labels_mapping <- importance_fid$tax
names(labels_mapping) <- importance_fid$fid
p2 <- ggplot(importance_fid, aes(x = importance_score, y = fid)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "CRC_16S", x = "Importance_score", y = "") +
  scale_y_discrete(labels = labels_mapping) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12))

### dietary fiber
importance_fid <- read.csv("../data/dietary_fiber/results/importance_fid.csv") %>% 
  arrange(importance_score) %>% tail(10)
colnames(importance_fid)[1] <- "fid"
importance_fid <- cbind(importance_fid, taxonomy_df_2[importance_fid$fid, ])
importance_fid$fid <- factor(importance_fid$fid, levels = importance_fid$fid)
importance_fid <- importance_fid %>%
  mutate(
    tax = case_when(
      Species != "s__" & !grepl("UMG|UBA|WWE", Species) ~ paste0(Species),
      Genus != "g__" & !grepl("UMG|UBA|WWE", Genus) ~ paste0(Genus),
      Family != "f__" & !grepl("UMG|UBA|WWE", Family) ~ paste0(Family),
      Order != "o__" & !grepl("UMG|UBA|WWE", Order) ~ paste0(Order),
      Class != "c__" & !grepl("UMG|UBA|WWE", Class) ~ paste0(Class),
      Phylum != "p__" & !grepl("UMG|UBA|WWE", Phylum) ~ paste0(Phylum),
      TRUE ~ "unclassified")
  )
labels_mapping <- ifelse(
  nchar(importance_fid$tax) > 20,
  paste0(substr(importance_fid$tax, 1, 20), "\n", substr(importance_fid$tax, 21, nchar(importance_fid$tax))),
  importance_fid$tax
)
names(labels_mapping) <- importance_fid$fid
p3 <- ggplot(importance_fid, aes(x = importance_score, y = fid)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Fiber_16S", x = "Importance_score", y = "") +
  scale_y_discrete(labels = labels_mapping) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12))

### CRC WGS
importance_fid <- read.csv("../data/CRC_WGS/results/importance_fid.csv") %>% 
  arrange(importance_score) %>% tail(10)
colnames(importance_fid)[1] <- "fid"
importance_fid$fid <- gsub("\\.", "_", importance_fid$fid)
importance_fid$fid <- ifelse(
  nchar(importance_fid$fid) > 20,
  paste0(substr(importance_fid$fid, 1, 20), "\n", substr(importance_fid$fid, 21, nchar(importance_fid$fid))),
  importance_fid$fid
)
importance_fid$fid <- factor(importance_fid$fid, levels = importance_fid$fid)
p4 <- ggplot(importance_fid, aes(x = importance_score, y = fid)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "CRC_WGS", x = "Importance_score", y = "") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12))

### T2D WGS
importance_fid <- read.csv("../data/T2D_WGS/results/importance_fid.csv") %>% 
  arrange(importance_score) %>% tail(10)
colnames(importance_fid)[1] <- "fid"
importance_fid$fid <- gsub("\\.", "_", importance_fid$fid)
importance_fid$fid <- ifelse(
  nchar(importance_fid$fid) > 20,
  paste0(substr(importance_fid$fid, 1, 20), "\n", substr(importance_fid$fid, 21, nchar(importance_fid$fid))),
  importance_fid$fid
)
importance_fid$fid <- factor(importance_fid$fid, levels = importance_fid$fid)
p5 <- ggplot(importance_fid, aes(x = importance_score, y = fid)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "T2D_WGS", x = "Importance_score", y = "") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12))

### Multi-class
importance_fid <- read.csv("../data/Multi_classification/results/importance_fid.csv") %>% 
  arrange(importance_score) %>% tail(10)
importance_fid$X <- as.character(importance_fid$X)
colnames(importance_fid)[1] <- "fid"
importance_fid <- cbind(importance_fid, taxonomy_df[importance_fid$fid, ])
importance_fid$fid <- factor(importance_fid$fid, levels = importance_fid$fid)
importance_fid <- importance_fid %>%
  mutate(
    tax = case_when(
      Species != "s__" ~ paste0(Species),
      Genus != "g__" ~ paste0(Genus),
      Family != "f__" ~ paste0(Family),
      Order != "o__" ~ paste0(Order),
      Class != "c__" ~ paste0(Class),
      Phylum != "p__" ~ paste0(Phylum),
      TRUE ~ "unclassified")
  )
labels_mapping <- importance_fid$tax
names(labels_mapping) <- importance_fid$fid
p6 <- ggplot(importance_fid, aes(x = importance_score, y = fid)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Multiclass", x = "Importance_score", y = "") +
  theme_bw() +
  scale_y_discrete(labels = labels_mapping) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 12),
        axis.text.y = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12))

p <- plot_grid(p1, p2, p3, p4, p5, align="hv", labels = c("a", "b", "c", "d", "e"),
               nrow = 2, ncol=3, plot=FALSE)
ggsave(p, filename = "Benchmark_tax.pdf",width = 14, height = 10, useDingbats=FALSE)
