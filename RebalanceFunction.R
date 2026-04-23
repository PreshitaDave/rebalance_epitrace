# ============================================================================
# EpiTrace Cell-Type Rebalancing Function (main function at line 48) 
# ============================================================================
.gini <- function(counts) {
  v <- sort(counts[counts > 0])
  n <- length(v)
  2 * sum(seq_len(n) * v) / (n * sum(v)) - (n + 1) / n
}

.print_summary <- function(alpha, mode, n_c, targets, classes, seed) {
  summary_df <- data.frame(
    celltype   = classes,
    n_before   = n_c[classes],
    n_after    = targets[classes],
    pct_change = round(100 * (targets[classes] / n_c[classes] - 1), 1),
    row.names  = NULL
  )
  summary_df <- summary_df[order(-summary_df$n_before), ]
  message(sprintf("\n── resample_cells  alpha=%.2f  mode='%s'  seed=%d ──",
                  alpha, mode, seed))
  message(sprintf("  Total cells : %d -> %d", sum(n_c), sum(targets)))
  message(sprintf("  Gini  before: %.3f   after: %.3f",
                  .gini(n_c), .gini(targets)))
  message(sprintf("  Ratio before: %.1fx  after: %.1fx",
                  max(n_c) / min(n_c), max(targets) / min(targets)))
  print(summary_df, row.names = FALSE)
}

.targets_down <- function(n_c, alpha) {
  n_min   <- min(n_c)
  targets <- round(n_c^(1 - alpha) * n_min^alpha)
  pmin(targets, n_c)
}

.targets_mixed <- function(n_c, alpha, cap_multiplier = 3) {
  G       <- exp(mean(log(n_c)))
  targets <- round(n_c^(1 - alpha) * G^alpha)
  cap     <- round(min(n_c) * cap_multiplier)
  ifelse(targets > n_c, pmin(targets, cap), targets)
}

.targets_up <- function(n_c, alpha) {
  n_max   <- max(n_c)
  targets <- round(n_c^(1 - alpha) * n_max^alpha)
  pmax(targets, n_c)
}

resample_cells <- function(seurat_obj,
                           alpha          = 0.5,
                           mode           = c("down", "mixed", "up"),
                           celltype_col   = "celltype",
                           cap_multiplier = 3,
                           seed           = 1234,
                           verbose        = TRUE) {
  
  mode <- match.arg(mode)
  stopifnot(
    is.numeric(alpha), length(alpha) == 1, alpha >= 0, alpha <= 1,
    celltype_col %in% colnames(seurat_obj@meta.data),
    is.numeric(cap_multiplier), cap_multiplier >= 1
  )
  
  set.seed(seed)
  
  meta <- seurat_obj@meta.data
  meta[[celltype_col]] <- as.character(meta[[celltype_col]])
  all_bcs   <- as.character(seurat_obj$cell)
  celltypes <- meta[[celltype_col]][match(all_bcs, as.character(meta$cell))]
  
  ct_tab     <- table(celltypes)
  classes    <- names(ct_tab)
  n_c        <- as.integer(ct_tab)
  names(n_c) <- classes
  
  targets <- switch(mode,
                    down  = .targets_down(n_c, alpha),
                    mixed = .targets_mixed(n_c, alpha, cap_multiplier),
                    up    = .targets_up(n_c, alpha)
  )
  names(targets) <- classes
  
  new_to_orig <- character(0)
  
  for (ct in classes) {
    pool <- all_bcs[celltypes == ct]
    tgt  <- targets[[ct]]
    
    if (tgt <= length(pool)) {
      drawn       <- sample(pool, size = tgt, replace = FALSE)
      new_to_orig <- c(new_to_orig, setNames(drawn, drawn))
      
    } else {
      gap   <- tgt - length(pool)
      extra <- sample(pool, size = gap, replace = TRUE)
      
      dup_counter <- list()
      extra_new   <- character(gap)
      extra_orig  <- character(gap)
      
      for (j in seq_len(gap)) {
        src <- extra[[j]]
        cnt <- dup_counter[[src]]
        cnt <- if (is.null(cnt) || is.na(cnt)) 1L else cnt + 1L
        dup_counter[[src]] <- cnt
        extra_new[[j]]  <- sprintf("%s_dup%d", src, cnt)
        extra_orig[[j]] <- src
      }
      
      new_to_orig <- c(new_to_orig,
                       setNames(pool,       pool),
                       setNames(extra_orig, extra_new))
    }
  }
  
  # ── Fast path: pure downsampling, no duplication ─────────────────────────
  needs_dup <- any(names(new_to_orig) != unname(new_to_orig))
  
  if (!needs_dup) {
    result_obj <- subset(seurat_obj, cells = names(new_to_orig))
    result_obj$original_cell  <- names(new_to_orig)[
      match(rownames(result_obj@meta.data), names(new_to_orig))]
    result_obj$resample_alpha <- alpha
    result_obj$resample_mode  <- mode
    if (verbose) .print_summary(alpha, mode, n_c, targets, classes, seed)
    return(result_obj)
  }
  
  # ── Slow path: upsampling required ───────────────────────────────────────
  orig_unique <- unique(unname(new_to_orig))
  base_obj    <- subset(seurat_obj, cells = orig_unique)
  
  orig_meta         <- base_obj@meta.data
  orig_meta_indexed <- orig_meta[match(unname(new_to_orig),
                                       as.character(orig_meta$cell)), , drop = FALSE]
  new_meta              <- orig_meta_indexed
  rownames(new_meta)    <- names(new_to_orig)
  new_meta$original_cell  <- unname(new_to_orig)
  new_meta$resample_alpha <- alpha
  new_meta$resample_mode  <- mode
  
  assay_names <- Seurat::Assays(base_obj)
  has_counts  <- vapply(assay_names, function(an) {
    mtx <- tryCatch(
      Seurat::GetAssayData(base_obj, assay = an, layer = "counts"),
      error = function(e) NULL)
    !is.null(mtx) && nrow(mtx) > 0 && ncol(mtx) > 0
  }, logical(1))
  assay_names_ok <- assay_names[has_counts]
  
  if (any(!has_counts))
    message("  Note: skipping assays with no counts layer: ",
            paste(assay_names[!has_counts], collapse = ", "))
  
  new_assay_list <- lapply(assay_names_ok, function(an) {
    mtx       <- Seurat::GetAssayData(base_obj, assay = an, layer = "counts")
    orig_cols <- match(unname(new_to_orig), colnames(mtx))
    new_mtx   <- mtx[, orig_cols, drop = FALSE]
    colnames(new_mtx) <- names(new_to_orig)
    Seurat::CreateAssayObject(counts = new_mtx,
                              min.cells = 0, min.features = 0,
                              check.matrix = FALSE)
  })
  names(new_assay_list) <- assay_names_ok
  
  # Start from base_obj skeleton and force-replace slots
  result_obj <- base_obj
  
  for (an in assay_names_ok) {
    result_obj@assays[[an]] <- new_assay_list[[an]]
  }
  
  # Restore assay keys (required by Seurat validator)
  for (an in assay_names_ok) {
    if (length(result_obj@assays[[an]]@key) == 0 ||
        result_obj@assays[[an]]@key == "") {
      result_obj@assays[[an]]@key <- paste0(tolower(gsub("[^a-zA-Z0-9]", "", an)), "_")
    }
  }
  
  # Set metadata with new unique barcodes as rownames
  result_obj@meta.data <- new_meta[names(new_to_orig), , drop = FALSE]
  
  # Reset active.ident to match new cell set
  result_obj@active.ident <- setNames(
    factor(rep(NA_character_, length(new_to_orig)),
           levels = levels(base_obj@active.ident)),
    names(new_to_orig)
  )
  
  # Clear stale slots
  result_obj@graphs     <- list()
  result_obj@neighbors  <- list()
  result_obj@reductions <- list()
  
  if (length(base_obj@reductions) > 0)
    message("  Note: reductions not transferred — ",
            "duplicated cells have no unique embedding.")
  
  if (verbose) .print_summary(alpha, mode, n_c, targets, classes, seed)
  return(result_obj)
}
