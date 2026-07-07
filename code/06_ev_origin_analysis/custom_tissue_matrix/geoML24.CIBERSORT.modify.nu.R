#' CIBERSORT R script v1.03
#' Note: Signature matrix construction is not currently available; use java version for full functionality.
#' Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' Requirements:
#'       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#'       install.packages('e1071')
#'       install.pacakges('parallel')
#'       install.packages('preprocessCore')
#'       if preprocessCore is not available in the repositories you have selected, run the following:
#'           source("http://bioconductor.org/biocLite.R")
#'           biocLite("preprocessCore")
#' Windows users using the R GUI may need to Run as Administrator to install or update packages.
#' This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
#' single-threaded in Windows.
#'
#' Usage:
#'       Navigate to directory containing R script
#'
#'   In R:
#'       source('CIBERSORT.R')
#'       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN)
#'
#'       Options:
#'       i)  perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#'       ii) QN = Quantile normalization of input mixture (default = TRUE)
#'
#' Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
#' Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
#' License: http://cibersort.stanford.edu/CIBERSORT_License.txt
#' Core algorithm
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @export
#' Modified CoreAlg function that returns results for all nu values
#' @param X cell-specific gene expression matrix
#' @param y mixed expression for a single sample
#' @param nu.v vector of nu values to try
#' @export
CoreAlg_multiNu <- function(X, y, nu.v = c(0.25, 0.5, 0.75)) {
  
  svn_itor <- length(nu.v)
  
  # Store results for each nu
  results_list <- list()
  nusvm <- rep(0, svn_itor)
  corrv <- rep(0, svn_itor)
  
  # Try each nu value
  for(t in 1:svn_itor) {
    nus <- nu.v[t]
    model <- svm(X, y, type = "nu-regression", kernel = "linear", nu = nus, scale = FALSE)
    
    # Calculate weights
    weights <- t(model$coefs) %*% model$SV
    weights[which(weights < 0)] <- 0
    w <- weights / sum(weights)
    
    # Reconstruct expression
    u <- sweep(X, MARGIN = 2, w, '*')
    k <- apply(u, 1, sum)
    
    # Calculate metrics
    nusvm[t] <- sqrt((mean((k - y)^2)))  # RMSE
    corrv[t] <- cor(k, y)  # Correlation
    
    # Store results for this nu
    results_list[[t]] <- list(
      model = model,
      w = w,
      rmse = nusvm[t],
      corr = corrv[t]
    )
  }
  
  # Find best model based on RMSE
  mn <- which.min(nusvm)
  best_model <- results_list[[mn]]$model
  best_w <- results_list[[mn]]$w
  best_rmse <- results_list[[mn]]$rmse
  best_corr <- results_list[[mn]]$corr
  
  # Return results for all nu values and the best one
  newList <- list(
    "all_results" = results_list,
    "best_w" = best_w,
    "best_rmse" = best_rmse,
    "best_corr" = best_corr,
    "best_nu" = nu.v[mn],
    "all_nu" = nu.v,
    "all_rmse" = nusvm,
    "all_corr" = corrv
  )
  
  return(newList)
}

#' Modified permutation function that selects best nu for each permutation
#' @param perm Number of permutations
#' @param X cell-specific gene expression
#' @param Y mixed expression matrix
#' @param nu.v vector of nu values to try
#' @export
doPerm_multiNu <- function(perm, X, Y, nu.v = c(0.25, 0.5, 0.75)) {
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  best_nu_dist <- rep(NA, perm)
  
  while(itor <= perm) {
    # Random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist), dim(X)[1])])
    
    # Standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    # Run modified core algorithm
    result <- CoreAlg_multiNu(X, yr, nu.v)
    
    # Store correlation from best nu
    mix_r <- result$best_corr
    best_nu_dist[itor] <- result$best_nu
    
    # Store correlation
    if(itor == 1) {
      dist <- mix_r
    } else {
      dist <- rbind(dist, mix_r)
    }
    
    itor <- itor + 1
  }
  
  newList <- list("dist" = dist, "best_nu_dist" = best_nu_dist)
  return(newList)
}

#' Main modified CIBERSORT function with per-sample nu optimization
#' @param sig_matrix file path to gene expression from isolated cells
#' @param mixture_file heterogenous mixed expression
#' @param perm Number of permutations (0 to skip)
#' @param QN Perform quantile normalization or not (TRUE/FALSE)
#' @param nu.v vector of nu values to try (default: c(0.25, 0.5, 0.75))
#' @export
CIBERSORT_perSampleNu <- function(sig_matrix, mixture_file, perm = 0, QN = TRUE, 
                                  nu.v = c(0.25, 0.5, 0.75)) {
  
  library(e1071)
  library(parallel)
  library(preprocessCore)
  
  # Read in data
  X <- read.table(sig_matrix, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  Y <- read.table(mixture_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  # Order by gene names
  X <- X[order(rownames(X)), ]
  Y <- Y[order(rownames(Y)), ]
  
  P <- perm  # Number of permutations
  
  # Anti-log if max < 50 in mixture file
  if(max(Y) < 50) {
    Y <- 2^Y
  }
  
  # Quantile normalization of mixture file
  if(QN == TRUE) {
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  # Intersect genes
  Xgns <- rownames(X)
  Ygns <- rownames(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX, ]
  XintY <- Xgns %in% rownames(Y)
  X <- X[XintY, ]
  
  # Standardize signature matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  # Empirical null distribution of correlation coefficients
  if(P > 0) {
    perm_results <- doPerm_multiNu(P, X, Y, nu.v)
    nulldist <- sort(perm_results$dist)
    best_nu_null <- perm_results$best_nu_dist
  } else {
    nulldist <- NULL
    best_nu_null <- NULL
  }
  
  # Prepare output
  header <- c('Mixture', colnames(X), 'Best_nu', 'P-value', 'Correlation', 'RMSE')
  
  # Additional columns for each nu value's metrics
  nu_headers <- c()
  for(nu in nu.v) {
    nu_headers <- c(nu_headers, paste0("RMSE_nu", nu), paste0("Corr_nu", nu))
  }
  
  # Update header
  header <- c(header[1:(length(header)-3)], nu_headers, 
              header[(length(header)-2):length(header)])
  
  output <- matrix()
  best_nu_per_sample <- rep(NA, ncol(Y))
  sample_names <- colnames(Y)
  
  # Iterate through mixtures
  for(itor in 1:ncol(Y)) {
    y <- Y[, itor]
    
    # Standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    # Run modified core algorithm
    result <- CoreAlg_multiNu(X, y, nu.v)
    
    # Get results
    w <- result$best_w
    mix_r <- result$best_corr
    mix_rmse <- result$best_rmse
    best_nu <- result$best_nu
    best_nu_per_sample[itor] <- best_nu
    
    # Collect metrics for all nu values
    all_nu_metrics <- c()
    for(i in 1:length(nu.v)) {
      all_nu_metrics <- c(all_nu_metrics, 
                          result$all_rmse[i], 
                          result$all_corr[i])
    }
    
    # Calculate p-value
    pval <- 9999
    if(P > 0 && !is.null(nulldist)) {
      pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))
    }
    
    # Prepare output row
    out <- c(sample_names[itor], w, best_nu, pval, mix_r, mix_rmse)
    out <- c(out[1:(length(out)-3)], all_nu_metrics, out[(length(out)-2):length(out)])
    
    if(itor == 1) {
      output <- out
    } else {
      output <- rbind(output, out)
    }
  }
  
  # Save results
  results_file <- "CIBERSORT-PerSampleNu-Results.txt"
  write.table(rbind(header, output), file = results_file, sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Prepare return object
  obj <- rbind(header, output)
  obj <- obj[-1, ]  # Remove header row
  obj_data <- matrix(as.numeric(unlist(obj[, -1])), nrow = nrow(obj))
  rownames(obj_data) <- obj[, 1]
  
  # Create column names
  cell_types <- colnames(X)
  metric_names <- c()
  for(nu in nu.v) {
    metric_names <- c(metric_names, paste0("RMSE_nu", nu), paste0("Corr_nu", nu))
  }
  
  colnames(obj_data) <- c(cell_types, "Best_nu", "P-value", "Correlation", "RMSE", metric_names)
  
  # Add additional information
  return_list <- list(
    results = obj_data,
    sample_best_nu = setNames(best_nu_per_sample, sample_names),
    nu_values = nu.v,
    perm_null_nu = best_nu_null
  )
  
  # If permutations were done, add null distribution
  if(P > 0) {
    return_list$null_distribution <- nulldist
  }
  
  # Also save summary statistics
  summary_stats <- data.frame(
    Sample = sample_names,
    Best_nu = best_nu_per_sample,
    Correlation = as.numeric(obj_data[, "Correlation"]),
    RMSE = as.numeric(obj_data[, "RMSE"]),
    P_value = as.numeric(obj_data[, "P-value"])
  )
  
  write.table(summary_stats, file = "CIBERSORT-PerSampleNu-Summary.txt", 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Print summary
  cat("CIBERSORT with per-sample nu optimization completed.\n")
  cat("Number of samples:", ncol(Y), "\n")
  cat("Nu values tested:", paste(nu.v, collapse = ", "), "\n")
  cat("Results saved to:", results_file, "\n")
  
  if(P > 0) {
    nu_table <- table(best_nu_per_sample)
    cat("\nBest nu distribution across samples:\n")
    for(nu in names(nu_table)) {
      cat("  nu =", nu, ":", nu_table[nu], "samples (", 
          round(100 * nu_table[nu] / length(sample_names), 1), "%)\n")
    }
    
    if(!is.null(best_nu_null)) {
      cat("\nBest nu distribution in null permutations:\n")
      null_nu_table <- table(best_nu_null)
      for(nu in names(null_nu_table)) {
        cat("  nu =", nu, ":", null_nu_table[nu], "permutations (", 
            round(100 * null_nu_table[nu] / P, 1), "%)\n")
      }
    }
  }
  
  return(return_list)
}

#' Wrapper function for backward compatibility
#' Uses the modified version with per-sample nu optimization
#' @export
CIBERSORT <- function(sig_matrix, mixture_file, perm = 0, QN = TRUE) {
  return(CIBERSORT_perSampleNu(sig_matrix, mixture_file, perm, QN))
}