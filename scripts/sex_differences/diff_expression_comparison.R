library(tibble)
library(dplyr)
library(DESeq2)
library(eulerr)
library(reshape2)
library(textshape)
library(ggplot2)

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/"
setwd(project_dir)

dds_kallisto <- readRDS("results/PRJNA285788/sex_differences/diff_expr/kallisto/dds_lrt.RDS")
dds_star <- readRDS("results/PRJNA285788/sex_differences/diff_expr/STAR/dds_lrt.RDS")

res_kallisto <- results(dds_kallisto)
res_star <- results(dds_star)

res_genes <- res_star %>% as_tibble() %>% mutate(GeneID = row.names(res_star)) %>% 
  dplyr::select(GeneID, baseMean:padj)

immune_genes <- read.table(
  "results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_genes.txt")[["V1"]]

canonical_immune_genes <- read.table(
  "results/PRJNA285788/immune_gene_identification/immune_genes/ob_bt_imm_genes.txt")[["V1"]]

# PCAs ----------------------------------------------------------------------------------------

pca_kallisto <- plotPCA(vst(dds_kallisto), intgroup = "Sex")
pca_star <- plotPCA(vst(dds_star), intgroup = "Sex")

ggsave("results/PRJNA285788/sex_differences/diff_expr/kallisto/pca.png", plot = pca_kallisto)
ggsave("results/PRJNA285788/sex_differences/diff_expr/STAR/pca.png", plot = pca_star)

# Overlap of DEGs between STAR and Kallisto ---------------------------------------------------

degs_overlap <- list()
degs_overlap[["kallisto"]] <- row.names(res_kallisto[!is.na(res_kallisto$padj) & 
                                                       res_kallisto$padj < 0.05, ])
degs_overlap[["star"]] <- row.names(res_star[!is.na(res_star$padj) & res_star$padj < 0.05, ])

deg_overlap_plot <- plot(eulerr::euler(degs_overlap), quantities = TRUE, 
                         main = "Overlap in DEGs (padj < 0.05)")

immune_degs_overlap <- list()
immune_degs_overlap[["kallisto"]] <- row.names(res_kallisto[!is.na(res_kallisto$padj) & 
                                                       res_kallisto$padj < 0.05 & 
                                                         row.names(res_kallisto) %in% immune_genes, ])
immune_degs_overlap[["star"]] <- row.names(res_star[!is.na(res_star$padj) & 
                                                      res_star$padj < 0.05 &
                                                      row.names(res_star) %in% immune_genes, ])

immune_deg_overlap_plot <- plot(eulerr::euler(immune_degs_overlap), quantities = TRUE, 
                         main = "Overlap in immune DEGs (padj < 0.05)")

ggsave(plot = deg_overlap_plot, filename = "results/PRJNA285788/sex_differences/diff_expr/deg_overlap.png", 
       device = "png")
ggsave(plot = immune_deg_overlap_plot, filename = "results/PRJNA285788/sex_differences/diff_expr/imm_deg_overlap.png", 
       device = "png")

# Immune DEG analysis -------------------------------------------------------------------------

res_degs <- res_genes %>% dplyr::filter(padj < 0.05)

res_immune_genes <- res_genes %>% dplyr::filter(GeneID %in% immune_genes)
res_immune_degs <- res_genes %>% dplyr::filter(GeneID %in% immune_genes, padj < 0.05)

vst_counts <- DESeq2::vst(DESeq2::counts(dds_star))

counts <- vst_counts[row.names(vst_counts) %in% immune_genes,
                                   c(paste0("SRR", 2895245:2895251))]

normalised_counts <- t(scale(t(counts), center = TRUE, scale = TRUE))
normalised_counts <- normalised_counts[complete.cases(normalised_counts), ]
normalised_counts_reordered <- normalised_counts[hclust(dist(normalised_counts))$order, ]

normalised_counts_tbl <- melt(normalised_counts_reordered) %>% 
  as_tibble() %>% 
  dplyr::rename("GeneID" = "Var1", "Sample" = "Var2", "ReadCount" = "value") %>% 
  dplyr::filter(GeneID %in% res_immune_degs$GeneID)

heatmap <- (ggplot(normalised_counts_tbl, aes(Sample, GeneID)) + 
  geom_tile(aes(fill = ReadCount), colour = "white") + 
  scale_fill_gradient(name = "Normalised \nRead Counts") + 
  scale_x_discrete(labels = c(paste0("Female ", 1:4), paste0("Male ", 1:3))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_blank(),
        aspect.ratio = 1, axis.title = element_blank(),
        axis.ticks.y = element_blank()))

ggsave("results/PRJNA285788/sex_differences/diff_expr/deg_heatmap.png", heatmap,
       width = 25, height = 15, unit = "cm")

res_immune_degs_canonical <- res_genes %>% dplyr::filter(GeneID %in% canonical_immune_genes, padj < 0.05)

vst_counts <- DESeq2::vst(DESeq2::counts(dds_star))

counts <- vst_counts[row.names(vst_counts) %in% canonical_immune_genes,
                     c(paste0("SRR", 2895245:2895251))]

normalised_counts_canonical <- t(scale(t(counts), center = TRUE, scale = TRUE))
normalised_counts_canonical <- normalised_counts[complete.cases(normalised_counts), ]
normalised_counts_reordered_canonical <- normalised_counts[hclust(dist(normalised_counts))$order, ]

normalised_counts_tbl_canonical <- melt(normalised_counts_reordered) %>% 
  as_tibble() %>% 
  dplyr::rename("GeneID" = "Var1", "Sample" = "Var2", "ReadCount" = "value") %>% 
  dplyr::filter(GeneID %in% res_immune_degs_canonical$GeneID)

heatmap_canonical <- (ggplot(normalised_counts_tbl_canonical, aes(Sample, GeneID)) + 
              geom_tile(aes(fill = ReadCount), colour = "white") + 
              scale_fill_gradient(name = "Normalised \nRead Counts") + 
              scale_x_discrete(labels = c(paste0("Female ", 1:4), paste0("Male ", 1:3))) +
              theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
                    axis.text.y = element_blank(),
                    aspect.ratio = 1, axis.title = element_blank(),
                    axis.ticks.y = element_blank()))

ggsave("results/PRJNA285788/sex_differences/diff_expr/deg_heatmap_canonical.png", heatmap_canonical,
       width = 25, height = 15, unit = "cm")


test_matrix <- matrix(data = c(length(res_immune_genes$GeneID), length(res_immune_degs$GeneID), 
                               length(res_genes$GeneID), length(res_degs$GeneID)), 
                             nrow = 2, ncol = 2)
f.test <- fisher.test(test_matrix)

out_file <- "results/PRJNA285788/sex_differences/diff_expr/stats.txt"
sink(out_file, append = FALSE)
cat(paste("Genes:", test_matrix[1,2], "\n"))
cat(paste("Immune Genes:", test_matrix[1,1], "\n"))
cat(paste("DEGs:", test_matrix[2,2], "\n"))
cat(paste("Immune DEGs:", test_matrix[2,1], "\n"))
cat("-------------------------------------------------", "\n")
cat(paste0("Share of immune genes diff. expr.: ", 
             round((test_matrix[2,1] / test_matrix[1,1]), 4) * 100, "%", "\n"))
cat(paste0("Share of all genes diff. expr.: ", 
             round((test_matrix[2,2] / test_matrix[1,2]), 4) * 100, "%", "\n"))
cat("-------------------------------------------------", "\n")
cat("\n")
cat(paste0("Fisher-test p-value: ", round(f.test$p.value, 3), "\n"))
sink()

save(heatmap, file = "results/PRJNA285788/sex_differences/diff_expr/heatmap.RData")


