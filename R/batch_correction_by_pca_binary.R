# Batch correction by PCA for RNA-seq DEG analysis: binary main variable
#
# Use this script when the main variable of interest is binary, for example:
#   - ALS case/control status
#   - exposure yes/no status
#
# Workflow:
#   1. Residualize normalized expression while protecting biological variables
#      of interest.
#   2. Extract RNA PCs from the residualized expression matrix.
#   3. Select RNA PCs carefully. Exclude PCs associated with the binary main
#      variable of interest using:
#        a) point-biserial correlation test, p < 0.05
#        b) Kruskal-Wallis test, p < 0.05
#        c) covariate-adjusted logistic regression, FDR < 0.05
#   4. Run limma-voom DEG analysis on the original count matrix, adding only
#      the retained RNA PCs as covariates.
#
# Expected data structure:
#   - counts: raw count matrix with genes in rows and samples in columns
#   - expr: normalized expression matrix, such as log2CPM, with genes in rows
#           and samples in columns
#   - meta: sample metadata with one row per sample

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
# Helper function: encode binary variable
# -----------------------------------------------------------------------------

encode_binary_outcome <- function(x, reference_level = NULL, case_level = NULL) {
  if (is.numeric(x)) {
    ux <- sort(unique(x[!is.na(x)]))
    if (!all(ux %in% c(0, 1))) {
      stop("Numeric binary variable must be coded as 0/1.")
    }
    return(as.numeric(x))
  }

  if (!is.null(reference_level) && !is.null(case_level)) {
    x_factor <- factor(x, levels = c(reference_level, case_level))
  } else {
    x_factor <- factor(x)
  }

  if (nlevels(x_factor) != 2) {
    stop(
      "Binary variable must have exactly two levels. Observed levels: ",
      paste(levels(x_factor), collapse = ", ")
    )
  }

  as.numeric(x_factor == levels(x_factor)[2])
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
# Diagnostic: scree plot to choose candidate number of PCs
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
# Step 2: select RNA PCs for a binary main variable of interest
# -----------------------------------------------------------------------------

run_pc_binary_association_tests <- function(meta,
                                            pc_cols,
                                            signal_col,
                                            reference_level = NULL,
                                            case_level = NULL) {
  missing_cols <- setdiff(c(pc_cols, signal_col), colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
  }

  signal_binary <- encode_binary_outcome(
    meta[[signal_col]],
    reference_level = reference_level,
    case_level = case_level
  )

  signal_factor <- factor(signal_binary, levels = c(0, 1))

  out <- lapply(pc_cols, function(pc) {
    pc_values <- meta[[pc]]

    point_biserial <- suppressWarnings(
      cor.test(signal_binary, pc_values, method = "pearson")
    )

    kruskal <- suppressWarnings(
      kruskal.test(pc_values ~ signal_factor)
    )

    data.frame(
      pc = pc,
      point_biserial_estimate = unname(point_biserial$estimate),
      point_biserial_p = point_biserial$p.value,
      kruskal_p = kruskal$p.value,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out)
}

run_pc_logistic_tests <- function(meta,
                                  pc_cols,
                                  outcome_col,
                                  covariates = c(),
                                  reference_level = NULL,
                                  case_level = NULL,
                                  adjust_method = "fdr") {
  missing_cols <- setdiff(c(pc_cols, outcome_col, covariates), colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
  }

  meta_model <- meta
  meta_model$.outcome <- encode_binary_outcome(
    meta[[outcome_col]],
    reference_level = reference_level,
    case_level = case_level
  )

  out <- lapply(pc_cols, function(pc) {
    model_terms <- c(pc, covariates)
    formula <- as.formula(paste(".outcome ~", paste(model_terms, collapse = " + ")))

    fit <- tryCatch(
      glm(formula, data = meta_model, family = binomial),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(data.frame(
        pc = pc,
        logistic_beta = NA_real_,
        logistic_p = NA_real_,
        model_status = "failed",
        stringsAsFactors = FALSE
      ))
    }

    coef_table <- summary(fit)$coefficients

    if (!pc %in% rownames(coef_table)) {
      return(data.frame(
        pc = pc,
        logistic_beta = NA_real_,
        logistic_p = NA_real_,
        model_status = "pc_dropped_from_model",
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      pc = pc,
      logistic_beta = coef_table[pc, "Estimate"],
      logistic_p = coef_table[pc, "Pr(>|z|)"],
      model_status = "ok",
      stringsAsFactors = FALSE
    )
  })

  result <- bind_rows(out)
  result$logistic_fdr <- p.adjust(result$logistic_p, method = adjust_method)
  result
}

select_rna_pcs_binary <- function(meta,
                                  pc_cols,
                                  signal_col,
                                  logistic_covariates = c(),
                                  reference_level = NULL,
                                  case_level = NULL,
                                  association_p_cutoff = 0.05,
                                  logistic_fdr_cutoff = 0.05,
                                  adjust_method = "fdr") {
  association_results <- run_pc_binary_association_tests(
    meta = meta,
    pc_cols = pc_cols,
    signal_col = signal_col,
    reference_level = reference_level,
    case_level = case_level
  )

  logistic_results <- run_pc_logistic_tests(
    meta = meta,
    pc_cols = pc_cols,
    outcome_col = signal_col,
    covariates = logistic_covariates,
    reference_level = reference_level,
    case_level = case_level,
    adjust_method = adjust_method
  )

  selection_table <- association_results %>%
    left_join(logistic_results, by = "pc") %>%
    mutate(
      exclude_by_point_biserial = point_biserial_p < association_p_cutoff,
      exclude_by_kruskal = kruskal_p < association_p_cutoff,
      exclude_by_logistic = logistic_fdr < logistic_fdr_cutoff,
      exclude_from_deg_model = exclude_by_point_biserial | exclude_by_kruskal | exclude_by_logistic,
      selection = ifelse(exclude_from_deg_model, "exclude", "include")
    )

  list(
    selection_table = selection_table,
    included_pcs = selection_table$pc[selection_table$selection == "include"],
    excluded_pcs = selection_table$pc[selection_table$selection == "exclude"],
    association_results = association_results,
    logistic_results = logistic_results
  )
}

# -----------------------------------------------------------------------------
# Step 3: limma-voom model with selected RNA PCs
# -----------------------------------------------------------------------------

run_limma_voom_with_rna_pcs <- function(counts,
                                        meta,
                                        sample_col = "Sample",
                                        group_col,
                                        reference_group,
                                        case_group,
                                        base_covariates = c(),
                                        rna_pcs = c(),
                                        cell_prop_covariates = c(),
                                        gene_annotation = NULL) {
  counts <- align_matrix_to_metadata(counts, meta, sample_col = sample_col)

  required_cols <- c(group_col, base_covariates, rna_pcs, cell_prop_covariates)
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
  }

  group <- factor(meta[[group_col]], levels = c(reference_group, case_group))
  if (any(is.na(group))) {
    stop("group_col contains values outside reference_group and case_group.")
  }

  y <- DGEList(
    counts = as.matrix(counts),
    group = group,
    genes = gene_annotation,
    samples = factor(meta[[sample_col]])
  )
  y <- calcNormFactors(y, method = "TMM")

  meta_model <- meta
  meta_model$group <- group

  model_terms <- c("group", base_covariates, rna_pcs, cell_prop_covariates)
  design_formula <- as.formula(paste("~", paste(model_terms, collapse = " + ")))
  design <- model.matrix(design_formula, data = meta_model)

  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  coef_name <- paste0("group", case_group)
  if (!coef_name %in% colnames(fit$coefficients)) {
    stop(
      "Could not find coefficient '", coef_name, "'. Available coefficients are: ",
      paste(colnames(fit$coefficients), collapse = ", ")
    )
  }

  deg <- topTable(fit, coef = coef_name, n = Inf, sort.by = "P")
  deg$diffexpressed <- "NO"
  deg$diffexpressed[deg$logFC > 0 & deg$adj.P.Val < 0.05] <- "UP"
  deg$diffexpressed[deg$logFC < 0 & deg$adj.P.Val < 0.05] <- "DOWN"

  list(
    design = design,
    voom = v,
    fit = fit,
    deg = deg,
    deg_summary = table(deg$diffexpressed)
  )
}

# -----------------------------------------------------------------------------
# Example usage for binary variable: exposure Yes/No or ALS case/control
# -----------------------------------------------------------------------------

# source("R/batch_correction_by_pca_binary.R")
#
# protected_covariates <- c(
#   "exposure_group",
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
# logistic_covariates <- c(
#   "age_at_sample_date", "sex", "subject_type", "Batch", "library_type",
#   "all_PC1", "all_PC2", "all_PC3", "all_PC4",
#   "CD4T_DNA", "CD8T_DNA", "Mono_DNA", "NK_DNA", "Bcell_DNA", "Neu_DNA"
# )
#
# pc_selection <- select_rna_pcs_binary(
#   meta = metadata_with_pcs,
#   pc_cols = candidate_pcs,
#   signal_col = "exposure_group",
#   logistic_covariates = logistic_covariates,
#   reference_level = "No",
#   case_level = "Yes",
#   association_p_cutoff = 0.05,
#   logistic_fdr_cutoff = 0.05,
#   adjust_method = "fdr"
# )
#
# selected_rna_pcs <- pc_selection$included_pcs
# write.csv(pc_selection$selection_table, "rna_pc_selection_binary.csv", row.names = FALSE)
#
# deg_results <- run_limma_voom_with_rna_pcs(
#   counts = count_matrix,
#   meta = metadata_with_pcs,
#   sample_col = "Sample",
#   group_col = "exposure_group",
#   reference_group = "No",
#   case_group = "Yes",
#   base_covariates = c(
#     "age_at_sample_date", "sex", "subject_type", "Batch", "library_type",
#     "all_PC1", "all_PC2", "all_PC3", "all_PC4"
#   ),
#   rna_pcs = selected_rna_pcs,
#   cell_prop_covariates = c(
#     "CD4T_DNA", "CD8T_DNA", "Mono_DNA", "NK_DNA", "Bcell_DNA", "Neu_DNA"
#   )
# )
