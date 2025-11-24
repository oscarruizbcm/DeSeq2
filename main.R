#!/usr/bin/env Rscript

suppressMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(matrixStats)
  library(tools)
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

# ------------------- Logging helper -------------------
log_file <- file.path(results_dir, "run_log.txt")

log_msg <- function(..., .type = "INFO") {
  ts  <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", .type, "] ", ts, " - ", paste(..., collapse = " "))
  # print to console (Docker logs)
  cat(msg, "\n")
  # append to log file
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("Starting DESeq2 pipeline.")
log_msg("INPUT_DIR =", counts_dir)
log_msg("OUTPUT_DIR =", results_dir)
log_msg("META_FILE =", meta_file)
log_msg("GROUP_COL =", group_col,
        "PADJ_CUTOFF =", padj_cutoff,
        "LFC_CUTOFF =", lfc_cutoff,
        "TOPN =", topN)

# ------------------- Wrap main logic in tryCatch -------------------
main <- function() {

  # ---------- Load metadata ----------
  log_msg("Step 1: Loading metadata from", meta_file)
  if (!file.exists(meta_file)) {
    stop("Metadata file does not exist: ", meta_file)
  }
  meta <- read.csv(meta_file, stringsAsFactors = FALSE)
  log_msg("Metadata loaded. Rows =", nrow(meta), "Cols =", ncol(meta))

  # ---------- Read & merge FeatureCounts outputs ----------
  log_msg("Step 2: Reading FeatureCounts files from", counts_dir)
  count_files <- list.files(counts_dir, pattern = "_counts\\.txt$", full.names = TRUE)
  if (length(count_files) == 0) stop("No _counts.txt files found in INPUT_DIR: ", counts_dir)
  log_msg("Found", length(count_files), "count files.")

  count_list <- lapply(count_files, function(f) {
    log_msg("Reading count file:", f)
    df <- read.table(f, header = TRUE, comment.char = "#", sep = "\t", check.names = FALSE)
    # Expect Geneid, Length, and a final column with counts
    if (!all(c("Geneid", "Length") %in% colnames(df))) {
      stop("File ", f, " is missing Geneid or Length columns.")
    }
    df <- df[, c("Geneid", "Length", colnames(df)[ncol(df)])]
    colnames(df)[3] <- file_path_sans_ext(basename(f))
    df
  })

  log_msg("Merging count files...")
  counts_merged <- Reduce(function(x, y) dplyr::full_join(x, y, by = c("Geneid", "Length")), count_list)
  counts_merged[is.na(counts_merged)] <- 0

  rownames(counts_merged) <- counts_merged$Geneid
  count_matrix <- counts_merged[, -(1:2), drop = FALSE]
  count_matrix <- as.matrix(count_matrix)
  mode(count_matrix) <- "numeric"
  count_matrix <- round(count_matrix)

  log_msg("Count matrix dimensions: genes =", nrow(count_matrix), "samples =", ncol(count_matrix))

  # ---------- Align metadata with count matrix ----------
  log_msg("Step 3: Aligning metadata with count matrix.")
  if (nrow(meta) != ncol(count_matrix)) {
    stop("Metadata row count (", nrow(meta),
         ") does not match number of samples (", ncol(count_matrix), ").")
  }
  rownames(meta) <- colnames(count_matrix)

  if (!(group_col %in% colnames(meta))) {
    stop("Metadata is missing the grouping column '", group_col, "'.")
  }

  meta[[group_col]] <- factor(meta[[group_col]])
  log_msg("Grouping column", group_col, "has levels:", paste(levels(meta[[group_col]]), collapse = ", "))

  # ---------- DESeq2 Analysis ----------
  log_msg("Step 4: Creating DESeqDataSet.")
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = meta,
                                design = as.formula(paste("~", group_col)))

  log_msg("Filtering low-count genes.")
  keep <- rowSums(counts(dds)) > 1
  dds <- dds[keep, ]
  log_msg("Remaining genes after filter:", nrow(dds))

  log_msg("Running DESeq().")
  dds <- DESeq(dds)
  log_msg("DESeq() complete.")

  # ---------- Normalized counts & global plots ----------
  log_msg("Step 5: Writing normalized counts.")
  norm_counts <- counts(dds, normalized = TRUE)
  write.csv(norm_counts, file = file.path(results_dir, "normalized_counts.csv"), quote = TRUE)

  log_msg("Step 6: PCA plot.")
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  pca_plot <- plotPCA(vsd, intgroup = group_col)
  pdf(file.path(results_dir, "PCA.pdf"), width = 7, height = 6)
  print(pca_plot)
  dev.off()

  log_msg("Step 7: Global heatmap of top variable genes.")
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

  # ---------- All pairwise contrasts ----------
  log_msg("Step 8: All pairwise contrasts.")
  levels_vec <- levels(meta[[group_col]])
  if (length(levels_vec) < 2) {
    stop("Grouping column '", group_col, "' has fewer than 2 levels.")
  }

  safe <- function(x) gsub("[^A-Za-z0-9_.-]+", "-", x)
  all_results <- list()

  for (i in 1:(length(levels_vec) - 1)) {
    for (j in (i + 1):length(levels_vec)) {
      A <- levels_vec[i]
      B <- levels_vec[j]
      label <- sprintf("%s_vs_%s", safe(A), safe(B))
      log_msg("Running contrast:", A, "vs", B, "(", label, ")")

      res <- results(dds, contrast = c(group_col, A, B))

      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      res_df$contrast <- paste(A, "vs", B)
      res_df$group_col <- group_col
      other_cols <- setdiff(colnames(res_df), c("contrast", "group_col", "gene"))
      res_df <- res_df[, c("contrast", "group_col", "gene", other_cols)]

      out_dir <- file.path(contrasts_dir, label)
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      write.csv(res_df, file = file.path(out_dir, paste0("DESeq2_", label, ".csv")), row.names = FALSE)

      # Significant genes for heatmap
      res_sig <- subset(res_df, !is.na(padj) & padj <= padj_cutoff & abs(log2FoldChange) >= lfc_cutoff)
      res_sig <- res_sig[order(res_sig$padj), , drop = FALSE]
      n_sig <- nrow(res_sig)
      log_msg("Contrast", label, "has", n_sig, "significant genes at padj <=", padj_cutoff,
              "and |LFC| >=", lfc_cutoff)

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
        note_file <- file.path(out_dir, paste0("NO_SIGNIFICANT_GENES_", label, ".txt"))
        writeLines(sprintf("No genes passed padj ≤ %.3f and |LFC| ≥ %.2f for contrast %s vs %s.",
                           padj_cutoff, lfc_cutoff, A, B),
                   con = note_file)
      }

      all_results[[label]] <- res_df
    }
  }

  log_msg("Step 9: Writing combined contrasts CSV.")
  combined <- do.call(rbind, all_results)
  write.csv(combined,
            file = file.path(results_dir, "DESeq2_all_pairwise_contrasts.csv"),
            row.names = FALSE)

  # Baseline contrast (first two levels)
  first_A <- levels_vec[1]; first_B <- levels_vec[2]
  log_msg("Step 10: Baseline contrast CSV for", first_A, "vs", first_B)
  res_first <- results(dds, contrast = c(group_col, first_A, first_B))
  res_first_df <- as.data.frame(res_first)
  res_first_df$gene <- rownames(res_first_df)
  res_first_df <- res_first_df[order(res_first_df$padj), ]
  write.csv(res_first_df,
            file = file.path(results_dir,
                             sprintf("DESeq2_results_%s_vs_%s.csv", safe(first_A), safe(first_B))),
            row.names = FALSE)

  log_msg("Pipeline finished successfully.")
}

# Actually run, with error logging
tryCatch(
  main(),
  error = function(e) {
    log_msg("ERROR:", conditionMessage(e), .type = "ERROR")
    # Capture and log a simple traceback if available
    tb <- paste(capture.output(traceback()), collapse = "\n")
    cat("[TRACEBACK]\n", tb, "\n", file = log_file, append = TRUE)
    stop(e)  # rethrow so Docker sees a failure exit code
  }
)
