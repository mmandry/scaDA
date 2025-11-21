sample.loc <- function(object,group.1,method=c("balanced","unbalanced")){
  set.seed(0810)
  # for one vs. all other case, use stratified sampling to get balanced sample size
  if(method=="balanced"){
    group.1.loc <- which(object@colData==group.1)
    num.group.1 <- length(group.1.loc)

    celltypes <- names(sort(table(object@colData), decreasing = TRUE))
    cellnums <- sort(table(object@colData), decreasing = TRUE)
    mat <- data.frame(celltype=celltypes,cellnum=as.numeric(cellnums))

    mat2 <- mat[mat$celltype!=group.1,]
    mat2$prop <- as.numeric(mat2$cellnum)/sum(as.numeric(mat2$cellnum))
    ncell <- nrow(mat2)

    group.2.loc <- NULL
    for (i in 1:ncell){
      cell <- mat2$celltype[i]
      cell.loc <- which(object@colData==cell)
      size <- floor(num.group.1 * mat2$prop[i])
      cell.loc.samp <- sample(cell.loc, size, replace = FALSE)
      group.2.loc <- c(group.2.loc, cell.loc.samp)
    }
  }

  if(method=="unbalanced"){
    group.1.loc <- which(object@colData == group.1)
    group.2.loc <- which(object@colData != group.1)
  }

  return(list(group.1.loc=group.1.loc, group.2.loc=group.2.loc))
}


#loglikelihood function for ZINB
zinb.loglink <- function(counts, p, u, k) {
  counts <- as.numeric(counts)
  n <- length(counts)
  
  # Pre-calculate common terms
  inv_k <- 1/k
  denom_term <- (1 / (1 + k * u))^(inv_k)
  
  # Initialize density vector
  dens <- numeric(n)
  
  # Identify zeros and non-zeros
  is_zero <- counts == 0
  
  # --- Logic for Zeros ---
  # P(Y=0) = p + (1-p) * (1 / (1+k*u))^(1/k)
  if (any(is_zero)) {
    dens[is_zero] <- p + (1 - p) * denom_term
  }
  
  # --- Logic for Non-Zeros ---
  # P(Y=y) = (1-p) * Gamma(...) * ...
  if (any(!is_zero)) {
    nz_counts <- counts[!is_zero]
    
    # Calculate Gamma term
    # g = gamma(y + 1/k) / (gamma(1/k) * gamma(y + 1))
    # Using lgamma is numerically safer and faster
    log_g <- lgamma(nz_counts + inv_k) - (lgamma(inv_k) + lgamma(nz_counts + 1))
    g <- exp(log_g)
    
    # Calculate probability
    term2 <- ((k * u) / (1 + k * u))^nz_counts
    dens[!is_zero] <- (1 - p) * g * denom_term * term2
  }
  
  # Sum of logs
  # Add small constant to avoid log(0) if density is effectively 0
  loglink <- sum(log(dens + 1e-300))
  
  return(loglink)
}

posterior.phi <- function(counts,p,u,k,m,tao){
  counts <- as.numeric(counts)
  zros <- counts[counts==0]
  nzros <- counts[counts!=0]
  n1 <- length(zros)
  n2 <- length(nzros)
  alpha <- 1/k
  temp1 <- 1/(1+k*u)
  get.phi <- n1*log(p+(1-p)*(temp1^alpha))+sum(log(gamma(nzros+alpha)))-n2*log(gamma(alpha))+
    n2*alpha*log(temp1)+sum(nzros*(log(k*u)-log(1+k*u)))-((log(k)-m)^2)/(2*(tao^2))-
    log(k)-log(tao)
  return(get.phi)
}

#### two kind of tasks: normal sample - human brain; disease sample - AD dataset

