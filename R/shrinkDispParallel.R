#' Shrink Dispersion Parameters (Parallelized & Robust)
#'
#' @param object scaDAdataset object
#' @param ncores Number of cores (default: detectCores() - 2)
#'
#' @return scaDAdataset object
#' @export
#' @import BiocParallel
#' @importFrom stats median IQR optimise
shrinkDispParallel <- function(object, ncores = NULL) {
  message("start shrink dispersion parameter")

  # --- 1. Setup Parallel Backend ---
  if (is.null(ncores)) {
    ncores <- parallel::detectCores() - 2
    if (is.na(ncores) || ncores < 1) ncores <- 1
  }
  
  # Use Chunking to prevent overhead
  n_tasks <- ncores * 50 
  message(paste("Using BiocParallel with", ncores, "workers and", n_tasks, "tasks."))
  
  BPPARAM <- BiocParallel::MulticoreParam(workers = ncores, tasks = n_tasks, progressbar = TRUE)

  count <- object@count
  peak_names <- rownames(count)
  group.1.loc <- object@params$g1
  group.2.loc <- object@params$g2
  
  # Keep sparse!
  dat <- count[, c(group.1.loc, group.2.loc)]
  nsam <- ncol(dat)
  nsam1 <- length(group.1.loc)
  
  poolCol <- 1:nsam
  cond1Col <- 1:nsam1
  cond2Col <- (nsam1 + 1):nsam
  npeak <- nrow(dat)

  # Helper to calculate empirical Bayes priors safely
  get_shrink_params <- function(phi_vec) {
    # Filter out NAs and infinite values
    valid_phi <- phi_vec[is.finite(phi_vec) & phi_vec > 0]
    
    if (length(valid_phi) < 10) {
      warning("Too few valid phi estimates to calculate shrinkage parameters.")
      return(list(m = 0, max.val = 100, tao = 1)) # Fallback defaults
    }
    
    log_inv_phi <- log(1/valid_phi)
    
    m <- median(log_inv_phi)
    
    # FIX: Handle max.val robustness
    # 1. Remove NAs
    # 2. Ensure it is at least slightly larger than the lower bound (0.01)
    raw_max <- max(1/valid_phi, na.rm = TRUE)
    max.val <- max(raw_max, 0.05) 
    
    sigma2.mar <- (IQR(log_inv_phi, na.rm=TRUE) / 1.349)^2
    tao <- sqrt(max(sigma2.mar * 0.5, 1e-2))
    
    return(list(m = m, max.val = max.val, tao = tao))
  }

  # Define Worker Function
  # Note: We pass the specific column indices and the specific param dataframe
  shrink_worker <- function(i, col_indices, param_df, prior_params) {
    counts <- as.numeric(dat[i, col_indices])
    
    # If previous step failed (NA), return NA
    if (is.na(param_df[i, ]$mu) || is.na(param_df[i, ]$p0)) {
      return(NA)
    }
    
    prev <- param_df[i, ]$p0
    nb_mu <- param_df[i, ]$mu
    
    # Use tryCatch to prevent one bad peak from crashing the job
    res <- tryCatch({
      # posterior.phi must be available in the environment or exported by scaDA
      # If it's internal, you might need scaDA:::posterior.phi
      optimise(posterior.phi, 
               c(0.001, prior_params$max.val), # Lowered bound slightly to 0.001
               tol = 1e-4, 
               maximum = TRUE, 
               counts = counts, 
               p = prev, 
               u = nb_mu, 
               m = prior_params$m, 
               tao = prior_params$tao)$maximum
    }, error = function(e) {
      return(param_df[i, ]$phi) # Return original (non-shrunk) phi on error
    })
    
    return(res)
  }

  # --- Run Phase 1: Pooled ---
  message("Shrinking phi for pooled data...")
  est_params_pooled <- object@params$param_pooled
  priors_pooled <- get_shrink_params(est_params_pooled$phi)
  
  res_pooled <- BiocParallel::bplapply(1:npeak, function(i) {
    shrink_worker(i, poolCol, est_params_pooled, priors_pooled)
  }, BPPARAM = BPPARAM)
  
  # Assign back (handling NAs if any remain)
  est_params_pooled$phi <- unlist(res_pooled)
  # Fill NAs with original values if shrinkage failed completely
  na_idx <- is.na(est_params_pooled$phi)
  if(any(na_idx)) est_params_pooled$phi[na_idx] <- object@params$param_pooled$phi[na_idx]


  # --- Run Phase 2: Condition 1 ---
  message("Shrinking phi for condition 1...")
  est_params_c1 <- object@params$param_g1
  priors_c1 <- get_shrink_params(est_params_c1$phi)
  
  res_c1 <- BiocParallel::bplapply(1:npeak, function(i) {
    shrink_worker(i, cond1Col, est_params_c1, priors_c1)
  }, BPPARAM = BPPARAM)
  
  est_params_c1$phi <- unlist(res_c1)
  na_idx <- is.na(est_params_c1$phi)
  if(any(na_idx)) est_params_c1$phi[na_idx] <- object@params$param_g1$phi[na_idx]


  # --- Run Phase 3: Condition 2 ---
  message("Shrinking phi for condition 2...")
  est_params_c2 <- object@params$param_g2
  priors_c2 <- get_shrink_params(est_params_c2$phi)
  
  res_c2 <- BiocParallel::bplapply(1:npeak, function(i) {
    shrink_worker(i, cond2Col, est_params_c2, priors_c2)
  }, BPPARAM = BPPARAM)
  
  est_params_c2$phi <- unlist(res_c2)
  na_idx <- is.na(est_params_c2$phi)
  if(any(na_idx)) est_params_c2$phi[na_idx] <- object@params$param_g2$phi[na_idx]


  # --- Update Object ---
  object@params <- list(g1 = group.1.loc,
                        g2 = group.2.loc,
                        param_pooled = est_params_pooled,
                        param_g1 = est_params_c1,
                        param_g2 = est_params_c2)
  
  message("Done.")
  return(object)
}
