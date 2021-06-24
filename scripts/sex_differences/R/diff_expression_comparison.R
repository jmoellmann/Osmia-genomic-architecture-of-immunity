library(DESeq2)
library(eulerr)
library(ggplot2)

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/"
setwd(project_dir)

dds_kallisto <- readRDS("results/PRJNA285788/sex_differences/diff_expr/kallisto/dds_lrt.RDS")
dds_star <- readRDS("results/PRJNA285788/sex_differences/diff_expr/STAR/dds_lrt.RDS")

res_kallisto <- results(dds_kallisto)
res_star <- results(dds_star)

write.csv(as.data.frame(res_kallisto), "results/PRJNA285788/sex_differences/diff_expr/kallisto/res.csv")
write.csv(as.data.frame(res_star), "results/PRJNA285788/sex_differences/diff_expr/STAR/res.csv")

pca_kallisto <- plotPCA(vst(dds_kallisto), intgroup = "Sex")
pca_star <- plotPCA(vst(dds_star), intgroup = "Sex")

ggsave("results/PRJNA285788/sex_differences/diff_expr/kallisto/pca.png", plot = pca_kallisto)
ggsave("results/PRJNA285788/sex_differences/diff_expr/STAR/pca.png", plot = pca_star)

degs <- list()
degs[["kallisto"]] <- row.names(res_kallisto[!is.na(res_kallisto$padj) & res_kallisto$padj < 0.05, ])
degs[["star"]] <- row.names(res_star[!is.na(res_star$padj) & res_star$padj < 0.05, ])

euler_plot <- plot(eulerr::euler(degs), quantities = TRUE, main = "Overlap in DEGs (padj < 0.05)")
ggsave(plot = euler_plot, filename = "results/PRJNA285788/sex_differences/diff_expr/euler_plot.png", 
       device = "png")
