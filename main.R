#!/usr/bin/env Rscript

suppressMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
})

# --- Environment variables ---
counts_dir <- Sys.getenv("INPUT_DIR", "/app/data/input")
results_dir <- Sys.getenv("OUTPUT_DIR", "/app/data/output")
meta_file <- Sys.getenv("META_FILE", "/app/data/metadata/sample_metadata.csv")

# Ensure results directory exists
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load metadata ---
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# --- Read and merge FeatureCounts outputs ---
library(dplyr)
library(tools)
count_files <- list.files(counts_dir, pattern = "_counts\\.txt$", full.names = TRUE)

count_list <- lapply(count_files, function(f) {
  df <- read.table(f, header = TRUE, comment.char = "#", sep = "\t", check.names = FALSE)
  df <- df[, c("Geneid", "Length", colnames(df)[ncol(df)])]
  colnames(df)[3] <- file_path_sans_ext(basename(f))
  df
})

counts_merged <- Reduce(function(x, y) full_join(x, y, by = c("Geneid", "Length")), count_list)

rownames(counts_merged) <- counts_merged$Geneid
count_matrix <- counts_merged[, -(1:2)]
count_matrix <- round(as.matrix(count_matrix))

# --- DESeq2 Analysis ---
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(matrixStats)   # for rowVars
})

# Make sure meta rows line up with columns of count_matrix
stopifnot(nrow(meta) == ncol(count_matrix))
rownames(meta) <- colnames(count_matrix)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = meta,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)

res_ordered <- res[order(res$padj), ]


# --- Save results to CSV ---
write.csv(as.data.frame(res_ordered),
          file = file.path(results_dir, "DESeq2_results.csv"),
          row.names = TRUE)

# --- Normalized counts ---
norm_counts <- counts(dds, normalized = TRUE)

# --- PCA on variance-stabilized data ---
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# plotPCA returns a ggplot; save to PDF
pca_plot <- plotPCA(vsd, intgroup = "condition")
pdf(file.path(results_dir, "PCA.pdf"), width = 7, height = 6)
print(pca_plot)
dev.off()

# --- Heatmap (top variable genes) ---
# Select top 50 most variable genes (adjust N as desired)
topN <- 50L
rv <- rowVars(norm_counts)
sel <- head(order(rv, decreasing = TRUE), topN)

# log2 transform for visualization
mat <- log2(norm_counts[sel, , drop = FALSE] + 1)

# Column annotations from metadata (expects a 'condition' column)
ann_col <- data.frame(condition = meta$condition)
rownames(ann_col) <- rownames(meta)

pdf(file.path(results_dir, "heatmap_top_variable_genes.pdf"), width = 8, height = 10)
pheatmap(mat,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = ann_col,
         show_rownames = FALSE,
         main = sprintf("Top %d most variable genes (log2 norm counts)", topN))
dev.off()


