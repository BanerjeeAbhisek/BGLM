# Graphical Lasso comparison for the prostate pathway analysis
#
# This script replaces the BGLM network-estimation step with graphical lasso
# (GLASSO). It uses the same prostate expression data and the same 16 KEGG
# pathways used in prostate.R.
#
# GLASSO estimates a sparse inverse covariance matrix for the gene-expression
# variables. It does not estimate the outcome-regression coefficients (Beta)
# produced by BVSGR(), and it does not use the PSA outcome Y. Therefore, this
# script is a direct comparator for the pathway dependence networks, not for
# the BGLM gene-outcome signal-screening or prediction results.
#
# Required packages:
# install.packages(c("glasso", "sna"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("KEGGgraph", "impute"))
#
# Run this file from the repository folder containing dat.csv, pathway.csv,
# myImagePlot.r, and extdata/.

required_packages <- c("glasso", "sna", "KEGGgraph", "impute")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

# gplot() looks up its layout functions on the attached package search path.
suppressPackageStartupMessages(library(sna))

required_paths <- c("dat.csv", "pathway.csv", "extdata")
missing_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_paths) > 0L) {
  stop(
    "Run this script from the Abhisek_BVSG repository folder. Missing: ",
    paste(missing_paths, collapse = ", ")
  )
}

# Prepare the same pathway-specific expression matrices as prostate.R

pathway <- read.csv(
  "pathway.csv",
  header = FALSE,
  stringsAsFactors = FALSE
)[[1]]

pathway <- pathway[-c(20:24, 29, 32, 34, 63, 64, 78, 170)]
xml_files <- file.path("extdata", paste0(pathway, ".xml"))

if (any(!file.exists(xml_files))) {
  stop("One or more KEGG XML files listed in pathway.csv are missing.")
}

pathway_gene_ids <- lapply(xml_files, function(xml_file) {
  graph <- KEGGgraph::parseKGML2Graph(xml_file, expandGenes = TRUE)
  KEGGgraph::nodes(graph)
})

pathway_sizes <- vapply(pathway_gene_ids, length, integer(1))
keep_pathways <- pathway_sizes < 1000L & pathway_sizes > 10L

subpathway <- pathway[keep_pathways]
subpathway_gene_ids <- pathway_gene_ids[keep_pathways]

# These are the same 16 pathway positions used in prostate.R.
ipath <- c(81, 82, 99, 103, 108, 111, 134, 136,
           145, 193, 195, 198, 199, 200, 205, 208)

pathway_names <- c(
  "MAPK Signaling Pathway",
  "ErbB Signaling Pathway",
  "mTOR Signaling Pathway",
  "Wnt Signaling Pathway",
  "Axon Guidance",
  "Focal Adhesion",
  "Long-term Potentiation",
  "Neurotrophin Signaling Pathway",
  "Insulin Signaling Pathway",
  "Pathways in Cancer",
  "Colorectal Cancer",
  "Endometrial Cancer",
  "Glioma",
  "Prostate Cancer",
  "Chronic Myeloid Leukemia",
  "Non-small Cell Lung Cancer"
)

if (max(ipath) > length(subpathway)) {
  stop("The pathway ordering no longer matches prostate.R.")
}

dat_raw <- read.csv(
  "dat.csv",
  header = FALSE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_ids <- paste0("hsa:", unlist(dat_raw[1, ], use.names = FALSE))
expression_raw <- as.matrix(dat_raw[-1, , drop = FALSE])
storage.mode(expression_raw) <- "numeric"
colnames(expression_raw) <- gene_ids

# impute.knn() expects genes in rows and subjects in columns.
expression_imputed <- t(impute::impute.knn(t(expression_raw))$data)
colnames(expression_imputed) <- gene_ids

# These are the subjects retained by the original prostate.R analysis.
expression_data <- expression_imputed[-c(15, 23, 47, 67, 82:101), , drop = FALSE]

make_pathway_matrix <- function(pathway_genes) {
  matched_columns <- match(pathway_genes, colnames(expression_data))
  matched_columns <- matched_columns[!is.na(matched_columns)]

  if (length(matched_columns) < 2L) {
    stop("A selected pathway has fewer than two genes in dat.csv.")
  }

  X <- expression_data[, matched_columns, drop = FALSE]

  # Remove any constant genes because a correlation matrix cannot be computed
  # for a zero-variance column.
  nonconstant <- apply(X, 2, stats::sd, na.rm = TRUE) > 0
  X <- X[, nonconstant, drop = FALSE]

  if (ncol(X) < 2L) {
    stop("A selected pathway has fewer than two nonconstant genes.")
  }

  X
}

selected_X <- lapply(subpathway_gene_ids[ipath], make_pathway_matrix)
names(selected_X) <- pathway_names
selected_codes <- subpathway[ipath]

# Fit GLASSO and select rho by five-fold cross-validation

fit_glasso_cv <- function(X, folds = 5L, number_rho = 40L,
                          edge_tolerance = 1e-6) {
  X_scaled <- scale(X)
  sample_covariance <- stats::cov(X_scaled)
  maximum_correlation <- max(abs(sample_covariance[upper.tri(sample_covariance)]))

  if (!is.finite(maximum_correlation) || maximum_correlation <= 0) {
    stop("Could not construct a valid rho sequence for a pathway.")
  }

  rho_grid <- exp(seq(
    log(max(maximum_correlation * 0.02, 1e-4)),
    log(maximum_correlation),
    length.out = number_rho
  ))

  # Resetting the seed gives every pathway the same patient split.
  set.seed(1)
  fold_id <- sample(rep(seq_len(folds), length.out = nrow(X)))
  validation_loss <- matrix(
    NA_real_,
    nrow = length(rho_grid),
    ncol = folds
  )

  for (fold in seq_len(folds)) {
    validation_rows <- fold_id == fold
    training_data <- X[!validation_rows, , drop = FALSE]
    validation_data <- X[validation_rows, , drop = FALSE]

    training_means <- colMeans(training_data)
    training_sds <- apply(training_data, 2, stats::sd)

    if (any(!is.finite(training_sds)) || any(training_sds == 0)) {
      stop("A gene is constant within a cross-validation fold.")
    }

    training_data <- scale(
      training_data,
      center = training_means,
      scale = training_sds
    )
    validation_data <- scale(
      validation_data,
      center = training_means,
      scale = training_sds
    )

    training_covariance <- stats::cov(training_data)
    validation_covariance <- stats::cov(validation_data)

    for (rho_index in seq_along(rho_grid)) {
      fit <- tryCatch(
        glasso::glasso(
          s = training_covariance,
          rho = rho_grid[rho_index],
          penalize.diagonal = FALSE,
          trace = FALSE
        ),
        error = function(e) NULL
      )

      if (is.null(fit) || fit$errflag != 0) {
        next
      }

      log_determinant <- determinant(fit$wi, logarithm = TRUE)

      if (log_determinant$sign <= 0) {
        next
      }

      # Negative Gaussian log-likelihood, excluding constants.
      validation_loss[rho_index, fold] <-
        sum(validation_covariance * fit$wi) -
        as.numeric(log_determinant$modulus)
    }
  }

  mean_loss <- rowMeans(validation_loss, na.rm = TRUE)

  if (!any(is.finite(mean_loss))) {
    stop("Cross-validation did not produce a valid GLASSO fit.")
  }

  valid_folds <- rowSums(is.finite(validation_loss))
  loss_se <- apply(validation_loss, 1, stats::sd, na.rm = TRUE) /
    sqrt(valid_folds)

  best_index <- which.min(mean_loss)
  selected_rho <- rho_grid[best_index]

  best_fit <- glasso::glasso(
    s = sample_covariance,
    rho = selected_rho,
    nobs = nrow(X),
    penalize.diagonal = FALSE,
    trace = FALSE
  )
  precision <- best_fit$wi

  # Convert the inverse covariance matrix to signed partial correlations.
  partial_correlation <- -precision /
    sqrt(outer(diag(precision), diag(precision)))
  diag(partial_correlation) <- 1
  dimnames(partial_correlation) <- list(colnames(X), colnames(X))
  dimnames(precision) <- list(colnames(X), colnames(X))

  adjacency <- abs(precision) > edge_tolerance
  diag(adjacency) <- FALSE
  dimnames(adjacency) <- dimnames(precision)

  list(
    rho = selected_rho,
    cv_loss = mean_loss[best_index],
    edge_count = sum(adjacency[upper.tri(adjacency)]),
    covariance = best_fit$w,
    precision = precision,
    partial_correlation = partial_correlation,
    adjacency = adjacency,
    cv_path = data.frame(
      rho = rho_grid,
      mean_validation_loss = mean_loss,
      standard_error = loss_se
    )
  )
}

glasso_results <- vector("list", length(selected_X))

for (pathway_index in seq_along(selected_X)) {
  message(
    sprintf(
      "Fitting GLASSO pathway %d of %d: %s",
      pathway_index,
      length(selected_X),
      pathway_names[pathway_index]
    )
  )

  glasso_results[[pathway_index]] <- fit_glasso_cv(
    selected_X[[pathway_index]],
    folds = 5L,
    number_rho = 40L
  )
}

names(glasso_results) <- pathway_names

# Save the numerical results and node-to-gene mapping

output_dir <- "glasso_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pathway_summary <- data.frame(
  figure = seq_along(pathway_names),
  pathway = pathway_names,
  KEGG_code = selected_codes,
  genes = vapply(selected_X, ncol, integer(1)),
  selected_rho = vapply(glasso_results, `[[`, numeric(1), "rho"),
  cross_validation_loss = vapply(glasso_results, `[[`, numeric(1), "cv_loss"),
  nonzero_edges = vapply(glasso_results, `[[`, integer(1), "edge_count"),
  plotted_edges = vapply(
    glasso_results,
    function(result) {
      sum(abs(result$partial_correlation[upper.tri(result$partial_correlation)]) >= 0.10)
    },
    integer(1)
  ),
  selection_rule = "5-fold CV, minimum validation loss"
)

utils::write.csv(
  pathway_summary,
  file.path(output_dir, "glasso_pathway_summary.csv"),
  row.names = FALSE
)

node_gene_mapping <- do.call(
  rbind,
  lapply(seq_along(selected_X), function(pathway_index) {
    data.frame(
      figure = pathway_index,
      pathway = pathway_names[pathway_index],
      KEGG_code = selected_codes[pathway_index],
      node = seq_len(ncol(selected_X[[pathway_index]])),
      gene = colnames(selected_X[[pathway_index]])
    )
  })
)

utils::write.csv(
  node_gene_mapping,
  file.path(output_dir, "glasso_node_gene_mapping.csv"),
  row.names = FALSE
)

saveRDS(glasso_results, file.path(output_dir, "glasso_results.rds"))

# Create one heatmap/network page for each of the 16 pathways

plot_glasso_heatmap <- function(partial_correlation, pathway_title) {
  matrix_to_plot <- abs(partial_correlation)
  p <- ncol(matrix_to_plot)
  label_positions <- unique(round(seq(1, p, length.out = min(p, 12L))))

  graphics::image(
    x = seq_len(p),
    y = seq_len(p),
    z = t(matrix_to_plot[p:1, , drop = FALSE]),
    col = gray.colors(10L, start = 1, end = 0.1),
    zlim = c(0, 1),
    axes = FALSE,
    xlab = "Gene node",
    ylab = "Gene node",
    main = paste0(pathway_title, "\nAbsolute GLASSO partial correlations")
  )

  graphics::axis(1, at = label_positions, labels = label_positions, cex.axis = 0.7)
  graphics::axis(
    2,
    at = p - label_positions + 1,
    labels = label_positions,
    las = 1,
    cex.axis = 0.7
  )
  graphics::box()
}

plot_glasso_network <- function(result, pathway_title, display_threshold = 0.10) {
  edge_weights <- abs(result$partial_correlation)
  edge_weights[edge_weights < display_threshold] <- 0
  diag(edge_weights) <- 0
  node_labels <- seq_len(nrow(edge_weights))
  displayed_edges <- sum(edge_weights[upper.tri(edge_weights)] > 0)

  if (displayed_edges == 0L) {
    graphics::plot.new()
    graphics::title(
      main = paste0(pathway_title, "\nGLASSO network"),
      sub = paste0("No edges with |partial correlation| >= ", display_threshold)
    )
    return(invisible(NULL))
  }

  sna::gplot(
    edge_weights,
    gmode = "graph",
    displaylabels = TRUE,
    label = node_labels,
    label.cex = 0.65,
    vertex.cex = 0.8,
    vertex.col = "red",
    edge.col = "blue",
    edge.lwd = 1 + 8 * edge_weights,
    thresh = 0,
    pad = 0.3,
    main = paste0(
      pathway_title,
      "\nGLASSO network",
      "\nrho = ",
      format(result$rho, digits = 3),
      ", edges shown = ",
      displayed_edges
    )
  )
}

grDevices::pdf(
  file.path(output_dir, "mygraph_glasso.pdf"),
  width = 12,
  height = 6
)
old_parameters <- graphics::par(no.readonly = TRUE)

for (pathway_index in seq_along(glasso_results)) {
  graphics::par(
    mfrow = c(1, 2),
    mar = c(4, 4, 5, 1),
    cex.main = 0.8
  )
  plot_glasso_heatmap(
    glasso_results[[pathway_index]]$partial_correlation,
    pathway_names[pathway_index]
  )
  plot_glasso_network(
    glasso_results[[pathway_index]],
    pathway_names[pathway_index]
  )
}

graphics::par(old_parameters)
grDevices::dev.off()

message("GLASSO analysis complete.")
message("Outputs saved in: ", normalizePath(output_dir))
message("Created glasso_output/mygraph_glasso.pdf")
message("Created glasso_output/glasso_pathway_summary.csv")
message("Created glasso_output/glasso_node_gene_mapping.csv")
message("Created glasso_output/glasso_results.rds")
