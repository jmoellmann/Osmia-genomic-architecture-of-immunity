library(biomaRt)
library(tximport)
library(dplyr)
library(stringr)
library(readr)
library(DESeq2)

de_analysis <- function(in_dir, metadata, tx2gene, design, test = "Wald", full = NULL, reduced = ~ 1){
  
  files <- file.path(in_dir, metadata$Run, "abundance.h5")
  
  txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
  txi_df <- data.frame(txi$counts)
  colnames(txi_df) <- metadata$Run
  
  dds <- DESeqDataSetFromTximport(txi, colData = metadata,
                                  design = design)
  
  if(test == "Wald"){
    return(DESeq(dds))
  } else if(test == "LRT"){
    return(DESeq(dds, test = test, full = full, reduced = reduced))
  }
}


# Start of Analysis ---------------------------------------------------------------------------

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics"

in_dir <- paste0(project_dir, "/results/PRJNA285788/sex_differences/quantification/kallisto")
out_dir <- paste0(project_dir, "/results/PRJNA285788/sex_differences/diff_expr/kallisto")

dir.create(out_dir, recursive = TRUE)


# get tx2gene information ---------------------------------------------------------------------


feature_table <- read.table(paste0(
  project_dir, "/data/reference/annotation/GCF_004153925.1_Obicornis_v3_feature_table.txt"), sep = "\t")

tx2gene <- feature_table[, c("V11", "V15")]

# get metadata --------------------------------------------------------------------------------

metadata <- read.csv(paste0(project_dir, "/data/datasets/PRJNA285788/accession/SraRunTable.csv"), 
                     header = TRUE)

experim_groups <- c("female", "male")

metadata <- metadata %>% dplyr::select(Run, LibrarySelection, sex) %>% 
  dplyr::rename(Sex = sex) %>% filter(LibrarySelection == "PolyA") %>% mutate_all(as.factor) 

# DEA, (LRT test) -------------------------------------------------------------------

dds_lrt <- de_analysis(in_dir, metadata, tx2gene, design = ~ Sex, test = "LRT", full = ~ Sex)


# save ----------------------------------------------------------------------------------------

saveRDS(dds_lrt, file = paste0(out_dir, "/dds_lrt.RDS"))