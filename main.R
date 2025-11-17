#!/usr/bin/env Rscript

suppressMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(matrixStats)
  library(tools)
  library(dplyr)
})

# ------------------- Environment variables -------------------
counts_dir  <- Sys.getenv("INPUT_DIR", "/app/data/input")
results_dir <- Sys.getenv("OUTPUT_DIR", "/app/data/output")
meta_file   <- Sys.getenv("META_FILE", "/app/data/metadata/sample_metadata.csv")

group_col   <- Sys.getenv("GROUP_COL", "condition")
padj_cutoff <- as.numeric(Sys.getenv("PADJ_CUTOFF", "0.05"))
lfc_cutoff  <- as.numeric(Sys.getenv("LFC_CUTOFF", "0"))
topN        <- as.integer(Sys.getenv("TOPN", "50"))

# Ensure directories
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
contrasts_dir <- file.path(results_dir, "contrasts")
if (!dir.exists(contrasts_dir)) dir.create(contrasts_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------- Load metadata -------------------
meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)

# ------------------- Read & merge FeatureCounts outputs -------------------
count_files <- list.files(counts_dir, pattern = "_counts\\.txt$", full.names = TRUE)
if (length(count_files) == 0) stop("No _counts.txt files found in INPUT_DIR.")

count_list <- lapply(count_files, function(f) {
  df <- read.table(f, header = TRUE, comment.char = "#", sep = "\t", check.names = FALSE)
  # Expect featureCounts columns include Geneid, Length, and a final column with counts
  df <- df[, c("Geneid", "Length", colnames(df)[ncol(df)])]
  colnames(df)[3] <- file_path_sans_ext(basename(f))
  df
})

counts_merged <- Reduce(function(x, y) full_join(x, y, by = c("Geneid", "Length")), count_list)
counts_merged[is.na(counts_merged)] <- 0

rownames(counts_merged) <- counts_merged$Geneid
count_matrix <- counts_merged[, -(1:2), drop = FALSE] %>% as.matrix()
mode(count_matrix) <- "numeric"
count_matrix <- round(count_matrix)

# ------------------- Align metadata with count matrix -------------------
# If metadata has a 'sample' column, use it to align; else assume row order matches columns.
if ("sample" %in% colnames(meta)) {
  # Force 'sample' to be character
  meta$sample <- as.character(meta$sample)
  # Check that all count columns have rows in meta
  missing_meta <- setdiff(colnames(count_matrix), meta$sample)
  if (length(missing_meta) > 0) {
    stop("These samples are in counts but missing from metadata$sample: ",
         paste(missing_meta, collapse = ", "))
  }
  meta <- meta %>% filter(sample %in% colnames(count_matrix)) %>%
    distinct(sample, .keep_all = TRUE) %>%
    column_to_rownames("sample")
  # Reorder metadata rows to match count columns
  meta <- meta[colnames(count_matrix), , drop = FALSE]
} else {
  if (nrow(meta) != ncol(count_matrix)) {
    stop("Metadata row count does not match number of samples, and no 'sample' column to align by.")
  }
  rownames(meta) <- colnames(count_matrix)
}

if (!(group_col %in% colnames(meta))) {
  stop(sprintf("Metadata is missing the grouping column '%s'.", group_col))
}

# Convert grouping column to factor with stable order
meta[[group_col]] <- factor(meta[[group_col]])

# ------------------- DESeq2 Analysis -------------------
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = meta,
                              design = as.formula(paste("~", group_col)))

# Filter out genes with almost all zeros (optional but common)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep, ]

dds <- DESeq(dds)

# Save normalized counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file = file.path(results_dir, "normalized_counts.csv"), quote = TRUE)

# Variance-stabilized data & PCA
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
pca_plot <- plotPCA(vsd, intgroup = group_col)
pdf(file.path(results_dir, "PCA.pdf"), width = 7, height = 6)
print(pca_plot)
dev.off()

# Heatmap of top variable genes (overall)
rv <- rowVars(norm_counts)
sel <- head(order(rv, decreasing = TRUE), topN)
mat_var <- log2(norm_counts[sel, , drop = FALSE] + 1)
ann_col <- meta[, group_col, drop = FALSE]
colnames(ann_col) <- group_col

pdf(file.path(results_dir, sprintf("heatmap_top_%d_variable_genes.pdf", topN)), width = 8, height = 10)
pheatmap(mat_var,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = ann_col,
         show_rownames = FALSE,
         main = sprintf("Top %d most variable genes (log2 norm counts)", topN))
dev.off()

# ------------------- All pairwise contrasts -------------------
levels_vec <- levels(meta[[group_col]])
if (length(levels_vec) < 2) {
  stop(sprintf("Grouping column '%s' has fewer than 2 levels.", group_col))
}

# Helper to sanitize file names
safe <- function(x) gsub("[^A-Za-z0-9_.-]+", "-", x)

# Collect all results across contrasts
all_results <- list()

for (i in 1:(length(levels_vec) - 1)) {
  for (j in (i + 1):length(levels_vec)) {
    A <- levels_vec[i]
    B <- levels_vec[j]
    label <- sprintf("%s_vs_%s", safe(A), safe(B))
    message("Running contrast: ", A, " vs ", B)

    res <- results(dds, contrast = c(group_col, A, B))
    # Add helpful columns
    res_tbl <- as.data.frame(res) %>%
      rownames_to_column("gene") %>%
      mutate(contrast = paste(A, "vs", B),
             group_col = group_col) %>%
      relocate(contrast, group_col, gene)

    # Save per-contrast CSV
    out_dir <- file.path(contrasts_dir, label)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(res_tbl, file = file.path(out_dir, paste0("DESeq2_", label, ".csv")), row.names = FALSE)

    # Per-contrast heatmap of top significant genes
    res_sig <- res_tbl %>% filter(!is.na(padj), padj <= padj_cutoff, abs(log2FoldChange) >= lfc_cutoff) %>%
      arrange(padj)
    n_sig <- nrow(res_sig)

    if (n_sig >= 2) {
      top_genes <- head(res_sig$gene, min(topN, n_sig))
      mat <- log2(norm_counts[top_genes, , drop = FALSE] + 1)

      pdf(file.path(out_dir, paste0("heatmap_", label, ".pdf")), width = 8, height = 10)
      pheatmap(mat,
               scale = "row",
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               annotation_col = ann_col,
               show_rownames = FALSE,
               main = sprintf("%s (top %d DE genes, padj ≤ %.3f, |LFC| ≥ %.2f)",
                              gsub("-", " ", label),
                              min(topN, n_sig), padj_cutoff, lfc_cutoff))
      dev.off()
    } else {
      # Write a small note if no significant genes
      writeLines(sprintf("No genes passed padj ≤ %.3f and |LFC| ≥ %.2f for contrast %s vs %s.",
                         padj_cutoff, lfc_cutoff, A, B),
                 con = file.path(out_dir, paste0("NO_SIGNIFICANT_GENES_", label, ".txt")))
    }

    all_results[[label]] <- res_tbl
  }
}

# Combined CSV of all contrasts
combined <- bind_rows(all_results, .id = "contrast_label")
write.csv(combined,
          file = file.path(results_dir, "DESeq2_all_pairwise_contrasts.csv"),
          row.names = FALSE)

# For convenience, also export a single "baseline" results table (first two levels)
first_A <- levels_vec[1]; first_B <- levels_vec[2]
res_first <- results(dds, contrast = c(group_col, first_A, first_B)) %>%
  as.data.frame() %>% rownames_to_column("gene")
write.csv(res_first[order(res_first$padj), ],
          file = file.path(results_dir, sprintf("DESeq2_results_%s_vs_%s.csv", safe(first_A), safe(first_B))),
          row.names = FALSE)

message("Analysis complete. Results in: ", results_dir)
