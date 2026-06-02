
# Save Illumina annotation for those CpGs that we have methylation data on
# Make gene tables: pull pathway annotation of relevance, make redundant where necessary (note unknown genes removed, less CpGs than in anno)
# This will then be used as the universe for permutation-based competitive testing

library(here)
library(tidyverse)
library(patchwork)
library(minfi)
library(IlluminaMouseMethylationanno.12.v1.mm10)
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("chiaraherzog/IlluminaMouseMethylationanno.12.v1.mm10")
library(stringr)
library(msigdbr)
library(biomaRt)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(officer)
library(flextable)

results_dir <- here("results")
dirOut <- file.path(results_dir,"0-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

# Illumina annotation CpGs #----
anno <- getAnnotation(IlluminaMouseMethylationanno.12.v1.mm10) %>% as.data.frame() #287,049 CpGs

anno <- anno %>%
  mutate(
    # Assign primary feature by priority: TSS200 > TSS1500 > Body
    primary_feature = case_when(
      grepl("tss_200", Feature_NCBI)    ~ "tss_200",
      grepl("tss_1500", Feature_NCBI)   ~ "tss_1500",
      grepl("tss_body", Feature_NCBI)   ~ "tss_body",
      TRUE                               ~ NA_character_
    ),
    
    # flag promoter CpGs (TSS200 or TSS1500)
    promoter = grepl("tss_200|tss_1500", Feature_NCBI)
  ) %>%
  mutate(
    gene_primary = sapply(strsplit(GeneName_NCBI, ";"), `[`, 1)
  )

load('data/beta_final.Rdata')

anno <- anno %>% dplyr::filter(Name %in% rownames(beta_final)) #285,698 CpGs
saveRDS(anno, file.path(dirOut, "anno_12.v1.mm10.Rds"))

rm(beta_final);gc()

# Gene tables #----

sum(is.na(anno$gene_primary)) #69,624
anno2 <- anno %>% filter(!is.na(gene_primary), gene_primary != "")

gene_tbl <- anno2 %>% distinct(gene_primary) #27,031 genes

# Standardize gene symbols #----
symbol_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(anno$gene_primary),
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
)

# keep only valid official symbols
symbol_map <- symbol_map %>%
  filter(!is.na(ENTREZID))

gene_tbl <- gene_tbl %>% left_join(symbol_map, by = c("gene_primary" = "SYMBOL"))
sum(is.na(gene_tbl$gene_primary)) # 0
anno2 <- anno2 %>% left_join(symbol_map, by = c("gene_primary" = "SYMBOL")) # 216,074 CpGs left
saveRDS(anno2, file = file.path(dirOut, "anno2_genes_filtered.Rds"))

# Annotate with msigdb #----
msig_mouse <- msigdbr(species = "Mus musculus")

# Hallmark, Oncogenic signatures, C2 Reactome
msig_mouse <- msig_mouse %>%
  filter(gs_collection %in% c("H", "C5", "C6") | gs_subcollection %in% c("CP:REACTOME")) %>%
  #dplyr::select(gs_collection, gs_subcollection, gs_collection_name, gs_name, gene_symbol, ensembl_gene)
  dplyr::select(gs_collection, gs_name, gene_symbol, ensembl_gene)

gene_tbl <- gene_tbl %>%
  left_join(msig_mouse,
            by = c("gene_primary" = "gene_symbol"))

# Hallmark - spermatogenesis #----
gene_H_tbl <- gene_tbl %>% dplyr::filter(gs_collection == "H") %>% dplyr::filter(gs_name != "HALLMARK_SPERMATOGENESIS")
table(gene_H_tbl$gs_name)
saveRDS(gene_H_tbl, file = file.path(dirOut, "hallmark_gene_universe.Rds"))

# Oncogenic signatures #----
gene_C6_tbl <- gene_tbl %>% dplyr::filter(gs_collection == "C6")
saveRDS(gene_C6_tbl, file = file.path(dirOut, "oncsign_gene_universe.Rds"))

# Reactome #----
gene_reactome_tbl <- gene_tbl %>% dplyr::filter(gs_collection == "C2")
saveRDS(gene_reactome_tbl, file = file.path(dirOut, "reactome_gene_universe.Rds"))
