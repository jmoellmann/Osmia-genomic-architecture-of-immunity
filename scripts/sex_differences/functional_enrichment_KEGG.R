library(fgsea)
library(DESeq2)
library(tidyverse)
library(KEGGREST)
library(EnrichmentBrowser)

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/"
setwd(project_dir)

ob2dm_gene <- read_csv("results/PRJNA285788/immune_gene_identification/ob_gene2dm_gene.csv")
  
dme_pathways <- getGenesets(org = "dme", db = "kegg", cache = TRUE, return.type="list")

res <- DESeq2::results(readRDS(paste0(
  "results/PRJNA285788/sex_differences/diff_expr/STAR/dds_lrt.RDS")))

res_tibble <- res %>% as_tibble() %>% mutate(Ob_GeneID = row.names(res)) %>% 
  left_join(ob2dm_gene) %>% mutate(Dm_GeneID = str_replace(Dm_GeneID, "[^0-9]+", "")) %>% 
  mutate(Dm_GeneID = str_replace(Dm_GeneID, "^0+", "")) %>% 
  filter(!is.na(Dm_GeneID))

stats <- res_tibble$pvalue
names(stats) <- res_tibble$Dm_GeneID
stats <- sort(stats)

fgsea_res <- fgsea(pathways = dme_pathways, stats = stats, scoreType = "std")
