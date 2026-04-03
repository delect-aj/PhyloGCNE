library(tidyverse)
library(reshape2)
library(cowplot)
library(patchwork)
library(stringr)
library(ggtext) 

setwd("/home/dongbiao/GCN/results")

### greengenes1 tax
tax_ref <- read.csv("/beegfs/db/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt", 
                    sep="\t", header = FALSE)
tax_ref$V1 <- as.character(tax_ref$V1)
taxonomy_df <- as.data.frame(str_split_fixed(tax_ref$V2, pattern = "; ", n = 7))
colnames(taxonomy_df) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(taxonomy_df) <- tax_ref$V1

### greengenes2 tax
tax_ref <- read.csv("/beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv", sep="\t")
taxonomy_df_2 <- as.data.frame(str_split_fixed(tax_ref$Taxon, pattern = "; ", n = 7))
colnames(taxonomy_df_2) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(taxonomy_df_2) <- tax_ref$Feature.ID

# Shared theme
theme_taxa <- function() {
  theme_bw() +
    theme(
      plot.title   = element_text(size = 14, face = "bold"),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_markdown(size = 14, hjust = 1),
      axis.title   = element_text(size = 14)
    )
}

# Label generator for GG1 (prefix format)
make_label_gg1 <- function(df) {
  df %>% mutate(
    tax = case_when(
      Species != "s__" ~ paste0("s__<i>", str_remove(Genus, "^g__"), " ",
                                str_remove(Species, "^s__"), "</i>"),
      Genus   != "g__" ~ paste0("g__<i>", str_remove(Genus,   "^g__"), "</i>"),
      Family  != "f__" ~ paste0("f__<i>", str_remove(Family,  "^f__"), "</i>"),
      Order   != "o__" ~ paste0("o__<i>", str_remove(Order,   "^o__"), "</i>"),
      Class   != "c__" ~ paste0("c__<i>", str_remove(Class,   "^c__"), "</i>"),
      Phylum  != "p__" ~ paste0("p__<i>", str_remove(Phylum,  "^p__"), "</i>"),
      TRUE ~ "unclassified"
    )
  )
}

# Label generator for GG2 (filter unknown taxa: UMG/UBA/WWE)
make_label_gg2 <- function(df) {
  df %>% mutate(
    tax = case_when(
      Species != "s__" & !grepl("UMG|UBA|WWE|GGB|SGB", Species) ~
        paste0("s__<i>", str_remove(Species, "^s__"), "</i>"),
      Genus   != "g__" & !grepl("UMG|UBA|WWE|GGB|SGB", Genus) ~
        paste0("g__<i>", str_remove(Genus,   "^g__"), "</i>"),
      Family  != "f__" & !grepl("UMG|UBA|WWE", Family) ~
        paste0("f__<i>", str_remove(Family,  "^f__"), "</i>"),
      Order   != "o__" & !grepl("UMG|UBA|WWE", Order) ~
        paste0("o__<i>", str_remove(Order,   "^o__"), "</i>"),
      Class   != "c__" & !grepl("UMG|UBA|WWE", Class) ~
        paste0("c__<i>", str_remove(Class,   "^c__"), "</i>"),
      Phylum  != "p__" & !grepl("UMG|UBA|WWE", Phylum) ~
        paste0("p__<i>", str_remove(Phylum,  "^p__"), "</i>"),
      TRUE ~ "unclassified"
    )
  )
}

## ── a: IBD 16S ──────────────────────────────────────────────────────────────
importance_fid <- read.csv("../data/IBD_16S/results/importance_fid.csv") %>%
  arrange(importance_score) %>% tail(10)
importance_fid$X <- as.character(importance_fid$X)
colnames(importance_fid)[1] <- "fid"
importance_fid <- cbind(importance_fid, taxonomy_df[importance_fid$fid, ])
importance_fid$fid <- factor(importance_fid$fid, levels = importance_fid$fid)
importance_fid <- make_label_gg1(importance_fid)
labels_mapping <- setNames(importance_fid$tax, importance_fid$fid)

p1 <- ggplot(importance_fid, aes(x = importance_score, y = fid)) +
  geom_segment(aes(x = 0, xend = importance_score, y = fid, yend = fid),
               color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue", fill = "white", shape = 21, stroke = 1) +
  labs(title = "IBD_16S", x = "Importance score", y = "") +
  scale_y_discrete(labels = labels_mapping) +
  theme_taxa()

## ── b: CRC 16S ──────────────────────────────────────────────────────────────
importance_fid <- read.csv("../data/CRC_16S/results/importance_fid.csv") %>%
  arrange(importance_score) %>% tail(10)
importance_fid$X <- as.character(importance_fid$X)
colnames(importance_fid)[1] <- "fid"
importance_fid <- cbind(importance_fid, taxonomy_df[importance_fid$fid, ])
importance_fid$fid <- factor(importance_fid$fid, levels = importance_fid$fid)
importance_fid <- make_label_gg1(importance_fid)
labels_mapping <- setNames(importance_fid$tax, importance_fid$fid)

p2 <- ggplot(importance_fid, aes(x = importance_score, y = fid)) +
  geom_segment(aes(x = 0, xend = importance_score, y = fid, yend = fid),
               color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue", fill = "white", shape = 21, stroke = 1) +
  labs(title = "CRC_16S", x = "Importance score", y = "") +
  scale_y_discrete(labels = labels_mapping) +
  theme_taxa()

## ── c: Fiber 16S ────────────────────────────────────────────────────────────
importance_fid <- read.csv("../data/dietary_fiber/results/importance_fid.csv") %>%
  arrange(importance_score) %>% tail(10)
colnames(importance_fid)[1] <- "fid"
importance_fid <- cbind(importance_fid, taxonomy_df_2[importance_fid$fid, ])
importance_fid$fid <- factor(importance_fid$fid, levels = importance_fid$fid)
labels_mapping <- setNames(
  c("s__<i>Prevotella disiens</i>", "g__<i>Bacteroides</i>", "g__<i>Ligilactobacillus</i>",
    "s__<i>Ruminococcus callidus</i>", "s__<i>Bacteroides uniformis</i>", "s__<i>Romboutsia ilealis</i>",
    "f__<i>Acutalibacteraceae</i>", "p__<i>Patescibacteria</i>", "s__<i>Akkermansia muciniphila</i>",
    "g__<i>Prevotella</i>"),
  importance_fid$fid
)

p3 <- ggplot(importance_fid, aes(x = importance_score, y = fid)) +
  geom_segment(aes(x = 0, xend = importance_score, y = fid, yend = fid),
               color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue", fill = "white", shape = 21, stroke = 1) +
  labs(title = "Fiber_16S", x = "Importance score", y = "") +
  scale_y_discrete(labels = labels_mapping) +
  theme_taxa()

## ── d: OSCC 16S ─────────────────────────────────────────────────────────────
importance_fid <- read.csv("../data/OSCC_16S/results/importance_fid.csv") %>%
  arrange(importance_score) %>% tail(10)
colnames(importance_fid)[1] <- "fid"
importance_fid <- cbind(importance_fid, taxonomy_df_2[importance_fid$fid, ])
importance_fid$fid <- factor(importance_fid$fid, levels = importance_fid$fid)
labels_mapping <- setNames(
  c("s__<i>Prevotella intermedia</i>", "s__<i>Prevotella intermedia</i>",
    "s__<i>Aggregatibacter aphrophilus</i>", "s__<i>Prevotella intermedia</i>",
    "s__<i>Aggregatibacter segnis</i>", "s__Neisseria sp000186165",
    "s__<i>Prevotella intermedia</i>", "g__<i>Haemophilus</i>",
    "s__<i>Neisseria meningitidis</i>", "s__<i>Neisseria perflava</i>"),
  importance_fid$fid
)

p4 <- ggplot(importance_fid, aes(x = importance_score, y = fid)) +
  geom_segment(aes(x = 0, xend = importance_score, y = fid, yend = fid),
               color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue", fill = "white", shape = 21, stroke = 1) +
  labs(title = "OSCC_16S", x = "Importance score", y = "") +
  scale_y_discrete(labels = labels_mapping) +
  theme_taxa()

## ── e: GC 16S ───────────────────────────────────────────────────────────────
importance_fid <- read.csv("../data/Gastritis_16S/results/importance_fid.csv") %>%
  arrange(importance_score) %>% tail(10)
colnames(importance_fid)[1] <- "fid"
importance_fid <- cbind(importance_fid, taxonomy_df_2[importance_fid$fid, ])
importance_fid$fid <- factor(importance_fid$fid, levels = importance_fid$fid)
labels_mapping <- setNames(
  c("s__<i>Gilliamella apicola</i> A", "s__Sutcliffiella_A horikoshii",
    "s__<i>Brachymonas chironomis</i>", "s__<i>Brachymonas denitrificans</i>",
    "s__<i>Anaerosinus glycerini</i>", "g__<i>Stenotrophomonas</i>",
    "s__<i>Brevundimonas aurantiaca</i>", "s__<i>Lactobacillus iners</i>",
    "s__<i>Lactobacillus hominis</i>", "g__<i>Lactobacillus</i>"),
  importance_fid$fid
)

p5 <- ggplot(importance_fid, aes(x = importance_score, y = fid)) +
  geom_segment(aes(x = 0, xend = importance_score, y = fid, yend = fid),
               color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue", fill = "white", shape = 21, stroke = 1) +
  labs(title = "GC_16S", x = "Importance score", y = "") +
  scale_y_discrete(labels = labels_mapping) +
  theme_taxa()

## ── f: CRC WGS ──────────────────────────────────────────────────────────────
importance_fid <- read.csv("../data/CRC_WGS/results/importance_fid.csv") %>%
  arrange(importance_score) %>% tail(10)
colnames(importance_fid)[1] <- "fid"
importance_fid$fid <- gsub("\\.", "_", importance_fid$fid)
importance_fid$fid <- factor(importance_fid$fid, levels = importance_fid$fid)
labels_mapping <- setNames(
  c("s__<i>Clostridiaceae</i> bacterium", "s__<i>Bacteroides clarus</i>",
    "s__<i>Ellagibacter isourolithinifaciens</i>", "s__<i>Clostridium</i> sp.",
    "s__<i>Prevotella copri</i>", "s__<i>Collinsella aerofaciens</i>",
    "s__<i>Phocaeicola plebeius</i>", "s__<i>Bacteroides cellulosilyticus</i>",
    "s__<i>Clostridiales</i> bacterium Choco116", "s__<i>Sutterella wadsworthensis</i>"),
  importance_fid$fid
)

p6 <- ggplot(importance_fid, aes(x = importance_score, y = fid)) +
  geom_segment(aes(x = 0, xend = importance_score, y = fid, yend = fid),
               color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue", fill = "white", shape = 21, stroke = 1) +
  labs(title = "CRC_WGS", x = "Importance score", y = "") +
  scale_y_discrete(labels = labels_mapping) +
  theme_taxa()

## ── g: T2D WGS ──────────────────────────────────────────────────────────────
importance_fid <- read.csv("../data/T2D_WGS/results/importance_fid.csv") %>%
  arrange(importance_score) %>% tail(10)
colnames(importance_fid)[1] <- "fid"
importance_fid$fid <- gsub("\\.", "_", importance_fid$fid)
importance_fid$fid <- factor(importance_fid$fid, levels = importance_fid$fid)
labels_mapping <- setNames(
  c("s__GGB6569 SGB9281", "s__<i>Bacteroides thetaiotaomicron</i>", "s__GGB9581 SGB4365",
    "s__<i>Lacrimispora celerecrescens</i>", "s__<i>Ruminococcus</i> sp. AF13-28",
    "s__<i>Rothia mucilaginosa</i>", "s__<i>Raoultibacter timonensis</i>",
    "s__Sanguibacter SGB15121", "s__<i>Lachnospiraceae</i> bacterium AM48-27BH",
    "s__<i>Dorea formicigenerans</i>"),
  importance_fid$fid
)

p7 <- ggplot(importance_fid, aes(x = importance_score, y = fid)) +
  geom_segment(aes(x = 0, xend = importance_score, y = fid, yend = fid),
               color = "steelblue", linewidth = 1) +
  geom_point(size = 3, color = "steelblue", fill = "white", shape = 21, stroke = 1) +
  labs(title = "T2D_WGS", x = "Importance score", y = "") +
  scale_y_discrete(labels = labels_mapping) +
  theme_taxa()

## ── Combined output ──────────────────────────────────────────────────────────
p <- plot_grid(p1, p2, p3, p4, p5, p6, p7,
               align = "hv", labels = c("a", "b", "c", "d", "e", "f", "g"),
               nrow = 3, ncol = 3, axis = "lb")

ggsave(p, filename = "../results/Benchmark_tax.pdf", width = 18, height = 12, useDingbats = FALSE)
ggsave(p, filename = "../results/Benchmark_tax.png", width = 18, height = 12)

