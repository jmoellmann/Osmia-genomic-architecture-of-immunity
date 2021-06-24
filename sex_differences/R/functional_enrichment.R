library(DESeq2)
library(ggplot2)
library(topGO)
library(dplyr)
library(stringr)
library(readr)
library(biomaRt)

ONTOLOGY = "CC"

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/"
setwd(project_dir)

source(paste0(project_dir, "scripts/general/R/GO_mappings_to_list.R"))

get_gsea_results <- function(stats, gene2GO, mode = "ks", alpha = 0.05, 
                             nodeSize = 20, topNodes = 200){
  
  if(mode == "fisher"){
    topGOobject <- new("topGOdata", ontology = ONTOLOGY,
                       allGenes = stats,
                       nodeSize = nodeSize, annot = annFUN.gene2GO, gene2GO = gene2GO)
    
    resultsFisher <- runTest(topGOobject, algorithm = "classic", statistic = "fisher")
    resultsFisher.weight01 <- runTest(topGOobject, algorithm = "weight01", statistic = "fisher")
    
    allRes <- GenTable(topGOobject, weight01Fisher = resultsFisher.weight01, 
                       classicFisher = resultsFisher, orderBy = "weight01Fisher", 
                       ranksOf = "weight01Fisher", topNodes = topNodes)
    allRes$weight01FisherPadj = p.adjust(allRes$weigh01Fisher, method = "BH")
    allRes$classicFisherPadj = p.adjust(allRes$classicFisher, method = "BH")
    
    
  } else if(mode == "ks"){
    topGOobject <- new("topGOdata", ontology = ONTOLOGY,
                       allGenes = stats, geneSel = function(x) x < alpha,
                       nodeSize = nodeSize, annot = annFUN.gene2GO, gene2GO = gene2GO)
    
    resultsKS <- runTest(topGOobject, algorithm = "classic", statistic = "ks")
    resultsKS.weight01 <- runTest(topGOobject, algorithm = "weight01", statistic = "ks")
    
    allRes <- GenTable(topGOobject, weight01KS = resultsKS.weight01, 
                       classicKS = resultsKS, orderBy = "weight01KS", 
                       ranksOf = "weight01KS", topNodes = topNodes)
    
    allRes$weight01KSpadj = p.adjust(allRes$weight01KS, method = "BH")
    allRes$classicKSpadj = p.adjust(allRes$classicKS, method = "BH")
  }
  
  return(allRes)
}


# Annotate O.bicornis genes with D.melanogaster GO terms --------------------------------------

ob_gene2ob_protein <- readr::read_csv("data/reference/annotation/geneID_to_proteinID.txt", 
                col_names = c("Ob_GeneID", "Ob_ProteinID")) %>% 
  mutate(Ob_GeneID = paste0("LOC", str_extract(Ob_GeneID, "[0-9]+")), 
         Ob_ProteinID = str_extract(Ob_ProteinID, "XP_[0-9]+.1"))

insect_orthologues_dir <- "results/PRJNA285788/phylogenetic_comparison/orthofinder_insects/Orthologues"
ob_protein2dm_protein <- readr::read_tsv(
  paste0(insect_orthologues_dir, 
         "/Orthologues_Osmia_bicornis/Osmia_bicornis__v__Drosophila_melanogaster.tsv")) %>% 
  tidyr::separate_rows(Osmia_bicornis, sep = ", ") %>% 
  tidyr::separate_rows(Drosophila_melanogaster, sep = ", ") %>% 
  dplyr::rename("Ob_ProteinID" = "Osmia_bicornis", "Dm_ProteinID" = "Drosophila_melanogaster") %>% 
  dplyr::select("Ob_ProteinID", "Dm_ProteinID")

bee_orthologues_dir <- "results/PRJNA285788/phylogenetic_comparison/orthofinder_bees/Orthologues"
ob_protein2bt_protein <- readr::read_tsv(
  paste0(bee_orthologues_dir, 
         "/Orthologues_Osmia_bicornis/Osmia_bicornis__v__Bombus_terrestris.tsv")) %>% 
  tidyr::separate_rows(Osmia_bicornis, sep = ", ") %>% 
  tidyr::separate_rows(Bombus_terrestris, sep = ", ") %>% 
  dplyr::rename("Ob_ProteinID" = "Osmia_bicornis", "Bt_ProteinID" = "Bombus_terrestris") %>% 
  dplyr::select("Ob_ProteinID", "Bt_ProteinID")

bt_mart <- biomaRt::useEnsemblGenomes(biomart = "metazoa_mart", dataset = "bterrestris_eg_gene")
bt_protein2bt_gene <- tibble::as_tibble(biomaRt::getBM(c("ensembl_peptide_id", "ensembl_gene_id"), 
                                                       mart = bt_mart)) %>% 
  dplyr::filter(ensembl_peptide_id != "") %>% 
  dplyr::rename("Bt_ProteinID" = "ensembl_peptide_id", "Bt_GeneID" = "ensembl_gene_id")

bt_protein2bt_GO <- biomaRt::getBM(attributes = c("ensembl_peptide_id", "go_id"), mart = bt_mart) %>% 
  as_tibble() %>% dplyr::rename("Bt_ProteinID" = "ensembl_peptide_id", "Bt_GO_ID" = "go_id") %>% 
  dplyr::filter("Bt_ProteinID" != "" & !is.na("Bt_ProteinID"))

# dm_mart <- biomaRt::useEnsemblGenomes(biomart = "metazoa_mart", dataset = "dmelanogaster_eg_gene")
# dm_protein2dm_GO <- biomaRt::getBM(attributes = c("refseq_peptide", "go_id"), mart = dm_mart) %>% 
#   as_tibble()

dm_protein2dm_GO <- readr::read_tsv("data/reference/annotation/D_melanogaster_proteinID_to_GO.tsv", 
                                 col_names = c("Dm_GO_ID", "Dm_ProteinID"), skip = 1) %>% 
  dplyr::select("Dm_ProteinID", "Dm_GO_ID")

ob_gene2dm_GO <- dplyr::left_join(ob_gene2ob_protein, dplyr::left_join(ob_protein2dm_protein, dm_protein2dm_GO)) %>% 
  filter(!is.na(Dm_GO_ID)) %>% dplyr::select("Ob_GeneID", "Dm_GO_ID")

ob_gene2dm_GO_list <- GO_mappings_to_list(as.data.frame(ob_gene2dm_GO))

ob_gene2bt_GO <- dplyr::left_join(ob_gene2ob_protein, dplyr::left_join(ob_protein2bt_protein, bt_protein2bt_GO)) %>% 
  filter(!is.na(Bt_GO_ID)) %>% dplyr::select("Ob_GeneID", "Bt_GO_ID")

ob_gene2bt_GO_list <- GO_mappings_to_list(as.data.frame(ob_gene2bt_GO))


bt_imm_orthl <- unique(as.vector(unlist(
  left_join(ob_gene2ob_protein, left_join(ob_protein2bt_protein, bt_protein2bt_gene)) %>% 
    filter(!is.na(Bt_GeneID), Bt_GeneID %in% read.table(
      "data/reference/putative_terrestris_immune_genes_ensembl.txt")[["V1"]]) %>% 
    dplyr::select(Ob_GeneID))))

dm_imm_orthl <- unique(as.vector(unlist(
  left_join(ob_gene2ob_protein, ob_protein2dm_protein) %>% 
    filter(!is.na(Dm_ProteinID), Dm_ProteinID %in% read.table(
      "data/reference/flybase_immune_process_list_proteins.txt")[["V1"]]) %>% 
    dplyr::select(Ob_GeneID))))

# Run GO enrichment analysis ------------------------------------------------------------------

dir.create(paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", ONTOLOGY), 
           showWarnings = FALSE, recursive = TRUE)
dir.create(paste0("results/PRJNA285788/sex_differences/diff_expr/kallisto/funct_enrich/", ONTOLOGY), 
           showWarnings = FALSE, recursive = TRUE)

for(tool in c("STAR", "kallisto")){
  res <- results(readRDS(paste0("results/PRJNA285788/sex_differences/diff_expr/", tool, "/dds_lrt.RDS")))
  
  stats <- res$padj
  names(stats) <- row.names(res)
  stats <- stats[!is.na(stats)]
  
  dm_imm_orthl_stats <- stats[names(stats) %in% dm_imm_orthl]
  bt_imm_orthl_stats <- stats[names(stats) %in% bt_imm_orthl]
    
  gsea_ns20 <- get_gsea_results(stats, ob_gene2dm_GO_list, nodeSize = 20, topNodes = 80)
  gsea_ns50 <- get_gsea_results(stats, ob_gene2dm_GO_list, nodeSize = 50, topNodes = 80)
  
  gsea_dm_imm_orthl_ns10 <- get_gsea_results(dm_imm_orthl_stats, ob_gene2dm_GO_list, nodeSize = 10,
                                             topNodes = 20)
  gsea_bt_imm_orthl_ns10 <- get_gsea_results(bt_imm_orthl_stats, ob_gene2dm_GO_list, nodeSize = 10,
                                             topNodes = 20)
  
  write_tsv(gsea_ns20, paste0("results/PRJNA285788/sex_differences/diff_expr/", 
                                   tool, "/funct_enrich/", ONTOLOGY, "/gsea_dm_ns20.tsv"))
  write_tsv(gsea_ns50, paste0("results/PRJNA285788/sex_differences/diff_expr/", 
                                   tool, "/funct_enrich/", ONTOLOGY, "/gsea_dm_ns50.tsv"))
  write_tsv(gsea_dm_imm_orthl_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/", 
                              tool, "/funct_enrich/", ONTOLOGY, "/gsea_dm_dm_imm_orthl_ns10.tsv"))
  write_tsv(gsea_bt_imm_orthl_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/", 
                              tool, "/funct_enrich/", ONTOLOGY, "/gsea_dm_bt_imm_orthl_ns10.tsv"))
  
}

for(tool in c("STAR", "kallisto")){
  res <- results(readRDS(paste0("results/PRJNA285788/sex_differences/diff_expr/", tool, "/dds_lrt.RDS")))
  
  stats <- res$padj
  names(stats) <- row.names(res)
  stats <- stats[!is.na(stats)]
  
  dm_imm_orthl_stats <- stats[names(stats) %in% dm_imm_orthl]
  bt_imm_orthl_stats <- stats[names(stats) %in% bt_imm_orthl]
  
  gsea_ns10 <- get_gsea_results(stats, ob_gene2bt_GO_list, nodeSize = 10, topNodes = 80)
  gsea_ns20 <- get_gsea_results(stats, ob_gene2bt_GO_list, nodeSize = 20, topNodes = 80)
  
  # adjust topNodes parameter down for ONTOLOGY = "CC" runs
  gsea_dm_imm_orthl_ns10 <- get_gsea_results(dm_imm_orthl_stats, ob_gene2bt_GO_list, nodeSize = 10, 
                                             topNodes = 15)
  gsea_bt_imm_orthl_ns10 <- get_gsea_results(bt_imm_orthl_stats, ob_gene2bt_GO_list, nodeSize = 10, 
                                             topNodes = 15)
  
  write_tsv(gsea_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/",
                                   tool, "/funct_enrich/", ONTOLOGY, "/gsea_bt_ns10.tsv"))
  write_tsv(gsea_ns20, paste0("results/PRJNA285788/sex_differences/diff_expr/", 
                             tool, "/funct_enrich/", ONTOLOGY, "/gsea_bt_ns20.tsv"))
  write_tsv(gsea_dm_imm_orthl_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/", 
                                           tool, "/funct_enrich/", ONTOLOGY, "/gsea_bt_dm_imm_orthl_ns10.tsv"))
  write_tsv(gsea_bt_imm_orthl_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/", 
                                           tool, "/funct_enrich/", ONTOLOGY, "/gsea_bt_bt_imm_orthl_ns10.tsv"))
  
}