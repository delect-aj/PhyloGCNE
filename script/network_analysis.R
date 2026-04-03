library(qiime2R)
library(ggClusterNet)
library(phyloseq)
library(igraph)
library(tidyverse)

setwd("/home/dongbiao/GCN/results")

metadata <- read.csv("../data/T2D_WGS/metadata.tsv", sep = "\t", row.names = 1) %>%
  mutate(Group = if_else(group == 0, "Control", "T2D"))
abundance_table <- read_q2biom("../data/T2D_WGS/table.biom") %>% t() %>%
  as.data.frame()
inter_sid <- intersect(rownames(abundance_table), rownames(metadata))
metadata <- metadata[inter_sid, ]
abundance_table <- abundance_table[inter_sid, ]

tax_df <- read.csv("/home/dongbiao/software/Phylo-Spec/database/WGS/wgs_mpa4_taxonomy.csv", 
                   check.names = FALSE) %>% distinct(Species, .keep_all = TRUE)
rownames(tax_df) <- gsub("_+", ".", tax_df$Species)
tax_df <- tax_df[colnames(abundance_table), ]

OTU <- otu_table(t(abundance_table), taxa_are_rows = TRUE)
SAM <- sample_data(metadata)
TAX <- tax_table(as.matrix(tax_df))
ps_raw <- phyloseq(OTU, SAM, TAX)

ps <- prune_taxa(taxa_sums(ps_raw) > 10, ps_raw)
# ps <- transform_sample_counts(ps_filtered, function(x) x / sum(x))

path_res <- "./network/"
dir.create(path_res, showWarnings = FALSE)
result <- network.2(
  ps           = ps,
  group        = "Group",
  fill         = "Phylum",
  N            = 0,                  
  big          = TRUE,
  maxnode      = 5,
  select_layout = TRUE,
  layout_net   = "model_maptree2",   
  r.threshold  = 0.5,
  p.threshold  = 0.05,
  label        = FALSE,
  path         = path_res,
  zipi         = FALSE,
  method       = "spearman"
)

p  <- result[[1]]       
tab <- result[[2]]

ggsave("../results/network.pdf", p, width = 8, height = 4,limitsize = FALSE)



