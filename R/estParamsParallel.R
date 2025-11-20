#' estimating ZINB parameters but using bplapply
#'
#' @param object scaDAdataset object
#' @param group.1 cell type 1
#' @param group.2 cell type 2
#' @param BPPARAM optional BiocParallel parameter (default = MulticoreParam())
#'
#' @return a scaDAdataset object contains parameter estimates in params slot
#' @export
#' @importFrom stats median plogis pchisq
#' @importFrom progress progress_bar
#' @importFrom BiocParallel bplapply MulticoreParam

estParamsParallel <- function(object, group.1=NULL, group.2=NULL, BPPARAM=BiocParallel::MulticoreParam(progressbar=TRUE)) {
  message("start initial parameter estimates")

  count <- object@count
  peak_names <- rownames(count)
  sfs <- apply(count, 2, mean)
  sfs <- sfs / median(sfs)

  if (is.null(group.2)) {
    cells_loc <- sample.loc(object=object, group.1=group.1, method="balanced")
    group.1.loc <- cells_loc$group.1.loc
    group.2.loc <- cells_loc$group.2.loc
  } else {
    group.1.loc <- which(object@colData == group.1)
    group.2.loc <- which(object@colData == group.2)
  }

  dat <- count[, c(group.1.loc, group.2.loc)]
  npeak <- nrow(dat)
  nsam1 <- length(group.1.loc)
  nsam2 <- length(group.2.loc)

  cond1Col <- seq_len(nsam1)
  cond2Col <- (nsam1 + 1):(nsam1 + nsam2)

  sfs_pooled <- sfs[c(group.1.loc, group.2.loc)]
  sfs_cell1 <- sfs[group.1.loc]
  sfs_cell2 <- sfs[group.2.loc]
  
  # Ensure zinb.loglink is available to workers
  # If it's internal to scaDA, you might need to export it or access via :::
  # zinb.loglink <- scaDA:::zinb.loglink 

  # Define worker function
  fit_one_peak <- function(i) {
    counts_pooled <- as.numeric(dat[i, ])
    counts_cell1 <- counts_pooled[cond1Col]
    counts_cell2 <- counts_pooled[cond2Col]

    ctrl <- pscl::zeroinfl.control(method = "L-BFGS-B")
    ctrl$reltol <- NULL
    ctrl$factr <- 1e-3/.Machine$double.eps

    tryCatch({
      # pooled
      m1 <- pscl::zeroinfl(counts ~ 1 + offset(log(sfs_pooled)) | 1 + offset(log(sfs_pooled)),
                           data = data.frame(counts = counts_pooled),
                           dist = "negbin", control = ctrl)
      mu1 <- exp(m1$coefficients$count)
      theta1 <- min(m1$theta, 150)
      p01 <- plogis(m1$coefficients$zero)
      
      # group1
      m2 <- pscl::zeroinfl(counts ~ 1 + offset(log(sfs_cell1)) | 1 + offset(log(sfs_cell1)),
                           data = data.frame(counts = counts_cell1),
                           dist = "negbin", control = ctrl)
      mu2 <- exp(m2$coefficients$count)
      theta2 <- min(m2$theta, 150)
      p02 <- plogis(m2$coefficients$zero)
      
      # group2
      m3 <- pscl::zeroinfl(counts ~ 1 + offset(log(sfs_cell2)) | 1 + offset(log(sfs_cell2)),
                           data = data.frame(counts = counts_cell2),
                           dist = "negbin", control = ctrl)
      mu3 <- exp(m3$coefficients$count)
      theta3 <- min(m3$theta, 150)
      p03 <- plogis(m3$coefficients$zero)

      # We return only what is needed for the params slot
      list(mu = c(mu1, mu2, mu3),
           theta = c(theta1, theta2, theta3),
           p0 = c(p01, p02, p03))
    }, error = function(e) {
      list(mu = rep(NA,3), theta = rep(NA,3), p0 = rep(NA,3))
    })
  }

  # Run parallel loop
  # Progress bar is handled by BPPARAM
  results <- BiocParallel::bplapply(seq_len(npeak), fit_one_peak, BPPARAM = BPPARAM)

  # Combine results
  est_params_pooled <- data.frame(mu = sapply(results, function(x) x$mu[1]),
                                  phi = sapply(results, function(x) x$theta[1]),
                                  p0 = sapply(results, function(x) x$p0[1]),
                                  row.names = peak_names)
  est_params_cell1 <- data.frame(mu = sapply(results, function(x) x$mu[2]),
                                 phi = sapply(results, function(x) x$theta[2]),
                                 p0 = sapply(results, function(x) x$p0[2]),
                                 row.names = peak_names)
  est_params_cell2 <- data.frame(mu = sapply(results, function(x) x$mu[3]),
                                 phi = sapply(results, function(x) x$theta[3]),
                                 p0 = sapply(results, function(x) x$p0[3]),
                                 row.names = peak_names)

  object@params <- list(g1 = group.1.loc,
                        g2 = group.2.loc,
                        param_pooled = est_params_pooled,
                        param_g1 = est_params_cell1,
                        param_g2 = est_params_cell2)

  return(object)
}
