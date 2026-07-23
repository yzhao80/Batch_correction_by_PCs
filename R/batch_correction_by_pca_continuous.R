# Batch correction by PCA for RNA-seq analysis: continuous main variable
#
# Use this script when the main variable of interest is continuous, for example:
#   - clinical score
#   - exposure burden score
#   - age-related phenotype
#
# Workflow:
#   1. Residualize normalized expression while protecting biological variables
#      of interest.
#   2. Extract RNA PCs from the residualized expression matrix.
#   3. Select RNA PCs carefully. Exclude PCs associated with the continuous main
#      variable of interest using:
#        a) Pearson correlation test, p < 0.05
#        b) Spearman correlation test, p < 0.05
#        c) covariate-adjusted linear regression, FDR < 0.05
#   4. Run a limma model on the original expression/count matrix, adding only
#      the retained RNA PCs as covariates.
#
# Note: logistic regression requires a binary outcome. For a continuous main
# variable, the covariate-adjusted analogue is linear regression. If the main
# variable is converted to a binary outcome, use batch_correction_by_pca_binary.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(limma)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Helper function: check and align expression/count matrix with metadata
# -----------------------------------------------------------------------------

align_matrix_to_metadata <- function(mat, meta, sample_col = "Sample") {
  if (!sample_col %in% colnames(meta)) {
    stop("sample_col was not found in metadata: ", sample_col)
  }

  missing_samples <- setdiff(meta[[sample_col]], colnames(mat))
  if (length(missing_samples) > 0) {
    stop(
      "The following metadata samples were not found in the matrix columns: ",
      paste(missing_samples, collapse = ", ")
    )
  }

  mat_aligned <- mat[, match(meta[[sample_col]], colnames(mat)), drop = FALSE]

  if (!all(meta[[sample_col]] == colnames(mat_aligned))) {
    stop("Metadata sample order does not match matrix column order after alignment.")
  }

  mat_aligned
}

# -----------------------------------------------------------------------------
# Step 1: residualize expression and extract RNA PCs
# -----------------------------------------------------------------------------

get_residual_expression_pcs <- function(expr,
                                        meta,
                                        sample_col = "Sample",
                                        protected_covariates,
                                        n_pcs = 20,
                                        center = TRUE,
                                        scale. = TRUE) {
  expr <- align_matrix_to_metadata(expr, meta, sample_col = sample_col)

  missing_covariates <- setdiff(protected_covariates, colnames(meta))
  if (length(missing_covariates) > 0) {
    stop(
      "The following protected covariates were not found in metadata: ",
      paste(missing_covariates, collapse = ", ")
    )
  }

  design_formula <- as.formula(
    paste("~ 1 +", paste(protected_covariates, collapse = " + "))
  )
  design <- model.matrix(design_formula, data = meta)

  # Each gene is modeled across samples. The resulting matrix has samples in
  # rows and genes in columns, which is the expected input format for prcomp().
  residual_exprs <- apply(expr, 1, function(gene_expr) {
    fit <- lm.fit(design, gene_expr)
    residuals(fit)
  })

  pca_result <- prcomp(residual_exprs, center = center, scale. = scale.)
  scores <- as.data.frame(pca_result$x)

  n_pcs <- min(n_pcs, ncol(scores))
  selected_scores <- scores[, seq_len(n_pcs), drop = FALSE]
  colnames(selected_scores) <- paste0("RNA_PC", seq_len(n_pcs))

  meta_with_pcs <- bind_cols(meta, selected_scores)

  list(
    residual_exprs = residual_exprs,
    pca_result = pca_result,
    scores = scores,
    meta_with_pcs = meta_with_pcs,
    pc_cols = colnames(selected_scores)
  )
}

# -----------------------------------------------------------------------------
# Optional diagnostic: scree plot to choose candidate number of PCs
# -----------------------------------------------------------------------------

plot_scree <- function(pca_result, max_pcs = 50) {
  pca_summary <- summary(pca_result)
  variance_explained <- pca_summary$importance[2, ]
  max_pcs <- min(max_pcs, length(variance_explained))

  scree_df <- data.frame(
    PC = seq_len(max_pcs),
    variance_explained = variance_explained[seq_len(max_pcs)]
  )

  ggplot(scree_df, aes(x = PC, y = variance_explained)) +
    geom_point() +
    geom_line() +
    theme_classic() +
    labs(
      x = "Principal component",
      y = "Proportion of variance explained",
      title = "PCA scree plot"
    )
}

# -----------------------------------------------------------------------------
# Step 2: select RNA PCs for a continuous main variable of interest
# -----------------------------------------------------------------------------

run_pc_continuous_correlation_tests <- function(meta,
                                                pc_cols,
                                                signal_col) {
  missing_cols <- setdiff(c(pc_cols, signal_col), colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
  }

  signal <- meta[[signal_col]]
  if (!is.numeric(signal)) {
    stop("For the continuous workflow, signal_col must be numeric.")
  }

  out <- lapply(pc_cols, function(pc) {
    pc_values <- meta[[pc]]

    pearson <- suppressWarnings(
      cor.test(signal, pc_values, method = "pearson")
    )

    spearman <- suppressWarnings(
      cor.test(signal, pc_values, method = "spearman", exact = FALSE)
    )

    data.frame(
      pc = pc,
      pearson_estimate = unname(pearson$estimate),
      pearson_p = pearson$p.value,
      spearman_estimate = unname(spearman$estimate),
      spearman_p = spearman$p.value,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out)
}

run_pc_linear_tests <- function(meta,
                                pc_cols,
                                outcome_col,
                                covariates = c(),
                                adjust_method = "fdr") {
  missing_cols <- setdiff(c(pc_cols, outcome_col, covariates), colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(meta[[outcome_col]])) {
    stop("For covariate-adjusted continuous testing, outcome_col must be numeric.")
  }

  out <- lapply(pc_cols, function(pc) {
    model_terms <- c(pc, covariates)
    formula <- as.formula(paste(outcome_col, "~", paste(model_terms, collapse = " + ")))

    fit <- tryCatch(
      lm(formula, data = meta),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(data.frame(
        pc = pc,
        linear_beta = NA_real_,
        linear_p = NA_real_,
        model_status = "failed",
        stringsAsFactors = FALSE
      ))
    }

    coef_table <- summary(fit)$coefficients

    if (!pc %in% rownames(coef_table)) {
      return(data.frame(
        pc = pc,
        linear_beta = NA_real_,
        linear_p = NA_real_,
        model_status = "pc_dropped_from_model",
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      pc = pc,
      linear_beta = coef_table[pc, "Estimate"],
      linear_p = coef_table[pc, "Pr(>|t|)"],
      model_status = "ok",
      stringsAsFactors = FALSE
    )
  })

  result <- bind_rows(out)
  result$linear_fdr <- p.adjust(result$linear_p, method = adjust_method)
  result
}

select_rna_pcs_continuous <- function(meta,
                                      pc_cols,
                                      signal_col,
                                      adjusted_covariates = c(),
                                      correlation_p_cutoff = 0.05,
                                      linear_fdr_cutoff = 0.05,
                                      adjust_method = "fdr") {
  correlation_results <- run_pc_continuous_correlation_tests(
    meta = meta,
    pc_cols = pc_cols,
    signal_col = signal_col
  )

  linear_results <- run_pc_linear_tests(
    meta = meta,
    pc_cols = pc_cols,
    outcome_col = signal_col,
    covariates = adjusted_covariates,
    adjust_method = adjust_method
  )

  selection_table <- correlation_results %>%
    left_join(linear_results, by = "pc") %>%
    mutate(
      exclude_by_pearson = pearson_p < correlation_p_cutoff,
      exclude_by_spearman = spearman_p < correlation_p_cutoff,
      exclude_by_adjusted_linear_model = linear_fdr < linear_fdr_cutoff,
      exclude_from_deg_model = exclude_by_pearson | exclude_by_spearman | exclude_by_adjusted_linear_model,
      selection = ifelse(exclude_from_deg_model, "exclude", "include")
    )

  list(
    selection_table = selection_table,
    included_pcs = selection_table$pc[selection_table$selection == "include"],
    excluded_pcs = selection_table$pc[selection_table$selection == "exclude"],
    correlation_results = correlation_results,
    linear_results = linear_results
  )
}

# -----------------------------------------------------------------------------
# Step 3a: limma-voom model with selected RNA PCs for count matrix input
# -----------------------------------------------------------------------------

run_limma_voom_continuous_with_rna_pcs <- function(counts,
                                                  meta,
                                                  sample_col = "Sample",
                                                  variable_col,
                                                  base_covariates = c(),
                                                  rna_pcs = c(),
                                                  cell_prop_covariates = c(),
                                                  gene_annotation = NULL) {
  counts <- align_matrix_to_metadata(counts, meta, sample_col = sample_col)

  required_cols <- c(variable_col, base_covariates, rna_pcs, cell_prop_covariates)
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
  }

  y <- DGEList(
    counts = as.matrix(counts),
    genes = gene_annotation,
    samples = factor(meta[[sample_col]])
  )
  y <- calcNormFactors(y, method = "TMM")

  model_terms <- c(variable_col, base_covariates, rna_pcs, cell_prop_covariates)
  design_formula <- as.formula(paste("~", paste(model_terms, collapse = " + ")))
  design <- model.matrix(design_formula, data = meta)

  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  if (!variable_col %in% colnames(fit$coefficients)) {
    stop(
      "Could not find coefficient '", variable_col, "'. Available coefficients are: ",
      paste(colnames(fit$coefficients), collapse = ", ")
    )
  }

  deg <- topTable(fit, coef = variable_col, n = Inf, sort.by = "P")
  deg$diffexpressed <- "NO"
  deg$diffexpressed[deg$logFC > 0 & deg$adj.P.Val < 0.05] <- "POSITIVE_ASSOCIATION"
  deg$diffexpressed[deg$logFC < 0 & deg$adj.P.Val < 0.05] <- "NEGATIVE_ASSOCIATION"

  list(
    design = design,
    voom = v,
    fit = fit,
    results = deg,
    result_summary = table(deg$diffexpressed)
  )
}

# -----------------------------------------------------------------------------
# Step 3b: limma model with selected RNA PCs for normalized expression input
# -----------------------------------------------------------------------------

run_limma_continuous_with_rna_pcs <- function(expr,
                                             meta,
                                             sample_col = "Sample",
                                             variable_col,
                                             base_covariates = c(),
                                             rna_pcs = c(),
                                             cell_prop_covariates = c()) {
  expr <- align_matrix_to_metadata(expr, meta, sample_col = sample_col)

  required_cols <- c(variable_col, base_covariates, rna_pcs, cell_prop_covariates)
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
  }

  model_terms <- c(variable_col, base_covariates, rna_pcs, cell_prop_covariates)
  design_formula <- as.formula(paste("~", paste(model_terms, collapse = " + ")))
  design <- model.matrix(design_formula, data = meta)

  fit <- lmFit(expr, design)
  fit <- eBayes(fit)

  if (!variable_col %in% colnames(fit$coefficients)) {
    stop(
      "Could not find coefficient '", variable_col, "'. Available coefficients are: ",
      paste(colnames(fit$coefficients), collapse = ", ")
    )
  }

  results <- topTable(fit, coef = variable_col, n = Inf, sort.by = "P")
  results$diffexpressed <- "NO"
  results$diffexpressed[results$logFC > 0 & results$adj.P.Val < 0.05] <- "POSITIVE_ASSOCIATION"
  results$diffexpressed[results$logFC < 0 & results$adj.P.Val < 0.05] <- "NEGATIVE_ASSOCIATION"

  list(
    design = design,
    fit = fit,
    results = results,
    result_summary = table(results$diffexpressed)
  )
}

# -----------------------------------------------------------------------------
# Example usage for continuous variable
# -----------------------------------------------------------------------------

# source("R/batch_correction_by_pca_continuous.R")
#
# protected_covariates <- c(
#   "clinical_score",
#   "sex",
#   "age_at_sample_date",
#   "subject_type"
# )
#
# pc_results <- get_residual_expression_pcs(
#   expr = log2cpm_matrix,
#   meta = metadata,
#   sample_col = "Sample",
#   protected_covariates = protected_covariates,
#   n_pcs = 20
# )
#
# metadata_with_pcs <- pc_results$meta_with_pcs
# plot_scree(pc_results$pca_result, max_pcs = 50)
# candidate_pcs <- paste0("RNA_PC", 1:12)
#
# adjusted_covariates <- c(
#   "age_at_sample_date", "sex", "subject_type", "Batch", "library_type",
#   "all_PC1", "all_PC2", "all_PC3", "all_PC4",
#   "CD4T_DNA", "CD8T_DNA", "Mono_DNA", "NK_DNA", "Bcell_DNA", "Neu_DNA"
# )
#
# pc_selection <- select_rna_pcs_continuous(
#   meta = metadata_with_pcs,
#   pc_cols = candidate_pcs,
#   signal_col = "clinical_score",
#   adjusted_covariates = adjusted_covariates,
#   correlation_p_cutoff = 0.05,
#   linear_fdr_cutoff = 0.05,
#   adjust_method = "fdr"
# )
#
# selected_rna_pcs <- pc_selection$included_pcs
# write.csv(pc_selection$selection_table, "rna_pc_selection_continuous.csv", row.names = FALSE)
#
# continuous_results <- run_limma_voom_continuous_with_rna_pcs(
#   counts = count_matrix,
#   meta = metadata_with_pcs,
#   sample_col = "Sample",
#   variable_col = "clinical_score",
#   base_covariates = c(
#     "age_at_sample_date", "sex", "subject_type", "Batch", "library_type",
#     "all_PC1", "all_PC2", "all_PC3", "all_PC4"
#   ),
#   rna_pcs = selected_rna_pcs,
#   cell_prop_covariates = c(
#     "CD4T_DNA", "CD8T_DNA", "Mono_DNA", "NK_DNA", "Bcell_DNA", "Neu_DNA"
#   )
# )
