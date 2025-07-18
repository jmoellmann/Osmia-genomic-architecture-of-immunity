library(DESeq2)
library(ggplot2)
library(topGO)
library(dplyr)
library(stringr)
library(readr)
library(biomaRt)

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/"
setwd(project_dir)

source("scripts/PRJNA285788/general/get_gsea_results.R")
# Read in immune gene orthologues -------------------------------------------------------------

immune_genes <- read.table(
  "results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_genes.txt")[["V1"]]

# Load Ob Gene to Dm/Bt GO terms tables -------------------------------------------------------

ob_gene2dm_GO <- read_csv("results/PRJNA285788/immune_gene_identification/funct_enrichment/ob_gene2dm_GO.csv")
ob_gene2bt_GO <- read_csv("results/PRJNA285788/immune_gene_identification/funct_enrichment/ob_gene2bt_GO.csv")

ob_gene2dm_GO_nona <- ob_gene2dm_GO %>% dplyr::filter(Ob_GeneID != "", Dm_GOID != "") %>% 
  distinct()
ob_gene2dm_GO_list <- GO_mappings_to_list(as.data.frame(ob_gene2dm_GO_nona))

ob_gene2bt_GO_nona <- ob_gene2bt_GO %>% dplyr::filter(Ob_GeneID != "", Bt_GOID != "") %>% 
  distinct()
ob_gene2bt_GO_list <- GO_mappings_to_list(as.data.frame(ob_gene2bt_GO_nona))


# Run GO enrichment analysis ------------------------------------------------------------------
for(treatment in c("imidacloprid", "thiacloprid")){
  for(ontology in c("CC", "MF", "BP")){
    
    out_dir <- paste0("results/PRJNA285788/pesticide_exposure/diff_splicing/", treatment, 
                      "/funct_enrich/", ontology)
    
    dir.create(out_dir, 
               showWarnings = FALSE, recursive = TRUE)
    
    res_all <- read_csv(paste0("results/PRJNA285788/pesticide_exposure/diff_splicing/", treatment, 
                        "/res_all.csv"))
    
    stats <- res_all$minPval
    names(stats) <- res_all$GeneID
    stats <- stats[!is.na(stats)]
    
    imm_gene_stats <- stats[names(stats) %in% immune_genes]
    
    gsea_ns20 <- get_gsea_results(stats, ob_gene2dm_GO_list, nodeSize = 20, ontology = ontology)
    gsea_ns50 <- get_gsea_results(stats, ob_gene2dm_GO_list, nodeSize = 50, ontology = ontology)
    gsea_ns100 <- get_gsea_results(stats, ob_gene2dm_GO_list, nodeSize = 100, ontology = ontology)
    gsea_ns150 <- get_gsea_results(stats, ob_gene2dm_GO_list, nodeSize = 150, ontology = ontology)
    
    gsea_bt_ns10 <- get_gsea_results(stats, ob_gene2bt_GO_list, nodeSize = 10, ontology = ontology)
    gsea_bt_ns20 <- get_gsea_results(stats, ob_gene2bt_GO_list, nodeSize = 20, ontology = ontology)
    
    gsea_imm_gene_ns10 <- get_gsea_results(imm_gene_stats, ob_gene2dm_GO_list, nodeSize = 10, 
                                           ontology = ontology)
    gsea_imm_gene_ns20 <- get_gsea_results(imm_gene_stats, ob_gene2dm_GO_list, nodeSize = 20, 
                                           ontology = ontology)
    gsea_imm_gene_ns50 <- get_gsea_results(imm_gene_stats, ob_gene2dm_GO_list, nodeSize = 50, 
                                           ontology = ontology)
    gsea_imm_gene_ns100 <- get_gsea_results(imm_gene_stats, ob_gene2dm_GO_list, nodeSize = 100, 
                                            ontology = ontology)
    gsea_imm_gene_ns150 <- get_gsea_results(imm_gene_stats, ob_gene2dm_GO_list, nodeSize = 150, 
                                            ontology = ontology)
    
    gsea_imm_gene_bt_ns10 <- get_gsea_results(imm_gene_stats, ob_gene2bt_GO_list, nodeSize = 10, 
                                              ontology = ontology)
    
    
    write_tsv(gsea_ns20, paste0(out_dir, "/gsea_dm_ns20.tsv"))
    write_tsv(gsea_ns50, paste0(out_dir, "/gsea_dm_ns50.tsv"))
    write_tsv(gsea_ns100, paste0(out_dir, "/gsea_dm_ns100.tsv"))
    write_tsv(gsea_ns150, paste0(out_dir, "/gsea_dm_ns150.tsv"))
    
    write_tsv(gsea_bt_ns10, paste0(out_dir, "/gsea_bt_ns10.tsv"))
    write_tsv(gsea_bt_ns20, paste0(out_dir, "/gsea_bt_ns20.tsv"))
    
    write_tsv(gsea_imm_gene_ns10, paste0(out_dir, "/gsea_imm_gene_ns10.tsv"))
    write_tsv(gsea_imm_gene_ns20, paste0(out_dir, "/gsea_imm_gene_ns20.tsv"))
    write_tsv(gsea_imm_gene_ns50, paste0(out_dir, "/gsea_imm_gene_ns50.tsv"))
    write_tsv(gsea_imm_gene_ns100, paste0(out_dir, "/gsea_imm_gene_ns100.tsv"))
    write_tsv(gsea_imm_gene_ns150, paste0(out_dir, "/gsea_imm_gene_ns150.tsv"))
    
    write_tsv(gsea_imm_gene_bt_ns10, paste0(out_dir, "/gsea_imm_gene_bt_ns10.tsv"))
    
  }  
}
