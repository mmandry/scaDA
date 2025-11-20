#' optimization of mu and phi parameters
#'
#' @param object scaDAdataset object
#' @param ncores Number of cores to use. Defaults to all cores - 2.
#'
#' @return scaDAdataset object
#' @export
#' @import BiocParallel
#' @importFrom stats optimise pchisq p.adjust
#'
optParamsParallel <- function(object, ncores = NULL) {
  message("start optimize parameter estimates")

  # --- 1. Setup Parallel Backend with CHUNKING ---
  if (is.null(ncores)) {
    ncores <- parallel::detectCores() - 2
    if (is.na(ncores) || ncores < 1) ncores <- 1 
  }
  
  # CRITICAL: Set tasks to reduce overhead. 
  # Splitting 200k peaks into ~2000 chunks prevents memory fragmentation.
  n_tasks <- ncores * 100
  
  message(paste("Using BiocParallel with", ncores, "workers and", n_tasks, "tasks."))
  BPPARAM <- BiocParallel::MulticoreParam(workers = ncores, 
                                          tasks = n_tasks,
                                          progressbar = TRUE)

  count <- object@count
  peak_names <- rownames(count)
  group.1.loc <- object@params$g1
  group.2.loc <- object@params$g2

  # Keep 'dat' SPARSE to save RAM (do not use as.matrix here for 80k cells)
  dat <- count[,c(group.1.loc, group.2.loc)]
  npeak <- dim(dat)[1]
  nsam <- dim(dat)[2]
  nsam1 <- length(group.1.loc)
  nsam2 <- length(group.2.loc)
  poolCol <- c(1:nsam)
  cond1Col <- c(1:nsam1)
  cond2Col <- c((nsam1+1):nsam)

  est_params_cell1 <- data.frame(object@params$param_g1)
  est_params_cell2 <- data.frame(object@params$param_g2)

  tol <- 1e-2 # Relaxed tolerance for speed
  nitr <- 5   # Reduced iterations for speed

  # Ensure zinb.loglink is available to workers
  # zinb.loglink <- scaDA:::zinb.loglink # Uncomment if package is not loaded

  # --- Helper Function for Optimization ---
  optimize_row <- function(i, col_indices, param_df) {
    # Explicitly convert sparse row to numeric vector to avoid S4 overhead in loop
    counts <- as.numeric(dat[i, col_indices])
    
    # Skip optimization if all counts are 0
    if(sum(counts) == 0) return(list(mu = 1e-6, prev = 1))

    prev <- param_df[i,]$p0
    nb_mu <- param_df[i,]$mu
    nb_phi <- param_df[i,]$phi
    
    # === FIX 1: CORRECT BOUNDS ===
    # Allow mu to go up to the max raw count, regardless of initialization
    current_max_count <- max(counts)
    max.mu <- max(max(param_df$mu), current_max_count + 10) 
    # =============================

    for (k in 1:nitr) {
      prev0 <- prev
      nb_mu0 <- nb_mu
      
      # Optimize Mu
      try({
        nb_mu <- optimise(zinb.loglink, c(0.01, max.mu), tol = 1e-4, maximum = TRUE, 
                          counts = counts, p = prev0, k = nb_phi)$maximum
      }, silent=TRUE)
      
      # Optimize Prev
      try({
        prev <- optimise(zinb.loglink, c(0.01, 1), tol = 1e-4, maximum = TRUE, 
                         counts = counts, u = nb_mu, k = nb_phi)$maximum
      }, silent=TRUE)

      if (abs(nb_mu0 - nb_mu) / (abs(nb_mu0)+1e-6) < tol && 
          abs(prev0 - prev) / (abs(prev0)+1e-6) < tol) {
        break
      }
    }
    return(list(mu = nb_mu, prev = prev))
  }

  # --- Run Optimization 1 ---
  message("Optimizing parameters for condition 1...")
  results_c1 <- BiocParallel::bplapply(1:npeak, function(i) {
    optimize_row(i, cond1Col, est_params_cell1)
  }, BPPARAM = BPPARAM)

  mu_opt_c1 <- sapply(results_c1, function(x) x$mu)
  prev_opt_c1 <- sapply(results_c1, function(x) x$prev)

  # --- Run Optimization 2 ---
  message("Optimizing parameters for condition 2...")
  results_c2 <- BiocParallel::bplapply(1:npeak, function(i) {
    optimize_row(i, cond2Col, est_params_cell2)
  }, BPPARAM = BPPARAM)

  mu_opt_c2 <- sapply(results_c2, function(x) x$mu)
  prev_opt_c2 <- sapply(results_c2, function(x) x$prev)

  # --- Statistics ---
  message("Calculating final statistics...")
  est_params_cell1$mu <- mu_opt_c1
  est_params_cell1$p0 <- prev_opt_c1
  est_params_cell2$mu <- mu_opt_c2
  est_params_cell2$p0 <- prev_opt_c2
  est_params_pooled <- object@params$param_pooled

  pval_zinb_shrink_opt <- numeric(npeak)
  tstats <- numeric(npeak)

  # Regular loop is fine for this part (fast vector math)
  for (i in 1:npeak){
    # Explicit numeric conversion again for safety
    counts_pool <- as.numeric(dat[i, poolCol])
    counts_c1   <- as.numeric(dat[i, cond1Col])
    counts_c2   <- as.numeric(dat[i, cond2Col])

    logL_null <- zinb.loglink(counts=counts_pool, p=est_params_pooled[i,]$p0, 
                              u=est_params_pooled[i,]$mu, k=est_params_pooled[i,]$phi)

    logL_alter_1 <- zinb.loglink(counts=counts_c1, p=est_params_cell1[i,]$p0, 
                                 u=est_params_cell1[i,]$mu, k=est_params_cell1[i,]$phi)

    logL_alter_2 <- zinb.loglink(counts=counts_c2, p=est_params_cell2[i,]$p0, 
                                 u=est_params_cell2[i,]$mu, k=est_params_cell2[i,]$phi)
    
    logL_alter <- logL_alter_1 + logL_alter_2

    stat_val <- -2*(logL_null - logL_alter)
    if(stat_val < 0) stat_val <- 0 # Prevent negative stats
    
    tstats[i] <- stat_val
    pval_zinb_shrink_opt[i] <- pchisq(stat_val, df=3, lower.tail = FALSE)
  }
  
  result <- data.frame(peakID=peak_names, tstats=tstats, pval=pval_zinb_shrink_opt)
  result$FDR <- p.adjust(result$pval, method='fdr')
  
  # Calculate fold change
  m1 <- (1 - est_params_cell1$p0) * est_params_cell1$mu
  m2 <- (1 - est_params_cell2$p0) * est_params_cell2$mu
  # Add small constant to avoid division by zero
  foch <- (m2 + 1e-6) / (m1 + 1e-6) 
  result$log2fc <- log2(foch)
  
  object@params <- list(g1 = group.1.loc,
                        g2 = group.2.loc,
                        param_pooled = object@params$param_pooled,
                        param_g1 = est_params_cell1,
                        param_g2 = est_params_cell2)
  object@result <- result
  
  message("Done.")
  return(object)
}
