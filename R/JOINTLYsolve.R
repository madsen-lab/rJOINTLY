#' @title JOINTsolve
#'
#' @description Function to solve the NMF problem in JOINTLY
#' 
#' @param kernel.list A list (per-batch) of kernels. See \link{alphaDecay}
#' @param snn.list A list (per-batch) of SNN graphs. See \link{computeSNN}
#' @param rare.list A list (per-batch) of rare cell indices. See \link{computeSNN}
#' @param cpca.result The results from CPCA analysis. See \link{cpca}
#' @param init The name of the method to use for initialization ("random" or "clustering"). [default = "clustering"]
#' @param norm.scale Boolean (TRUE or FALSE) indicating whether or not to scale the normalized counts [default = TRUE]
#' @param norm.minmax Boolean (TRUE or FALSE) indicating whether or not to MinMax scale the normalized counts [default = FALSE]
#' @param norm.center Boolean (TRUE or FALSE) indicating whether or not to center the normalized counts [default = FALSE]
#' @param k The number of factors to calculate [default = 15]
#' @param m The fuzzy parameter of for cluster initialization [default = 2]
#' @param iter.max The maximum number of iterations to perform [default = 200]
#' @param alpha Alpha parameter of the loss function [default = 1]
#' @param mu Mu parameter of the loss function [default = 20]
#' @param lambda Lambda parameter of the loss function [default = 1]
#' @param beta Beta parameter of the loss function [default = 5]
#' @param progressbar A logical indicating whether or not to print a progress bar [default = TRUE]
#' @param share.objects Boolean (TRUE or FALSE) indicating to use shared object to reduce memory requirements for parallel processing [default = TRUE]
#' @param ncpu The number of cpus to use for matrix multiplication [default = 1]
#' @param save_all Boolean (TRUE or FALSE) indicating if all H, F and W matrices should be saved and outputted [default = FALSE]
#'
#' @return A list containing normalized H matrices to use for clustering, as well a list (per-batch) of H matrices, a list (per-batch) of F matrices and a list (per-batch) of W matrices.
#' @import future
#' @import future.apply
#' @import e1071
#' @import SharedObject
#' @export

JOINTLYsolve <- function(kernel.list, snn.list, rare.list, cpca.result, init = "clustering", norm.scale = TRUE, norm.minmax = FALSE, norm.center = FALSE, k = 15, m = 2, iter.max = 200, alpha = 1, mu = 20, lambda = 1, beta = 5, progressbar = TRUE, share.objects = TRUE, ncpu = 1, save_all = FALSE) {
  ## Convert to dense matrices
  norm.list <- list()
  for (ds in 1:length(kernel.list)) { 
    kernel.list[[ds]] <- as.matrix(kernel.list[[ds]]) 
    snn.list[[ds]] <- as.matrix(snn.list[[ds]]) 
    norm.list[[ds]] <- as.matrix(cpca.result$normalized[[ds]]) 
  }
  
  ## Scale the normalized matrix if requested
  if (norm.scale) {
    for (ds in 1:length(kernel.list)) { 
      norm.list[[ds]] <- t(scale(t(norm.list[[ds]]), center = norm.center))
      if (norm.center) {
        norm.list[[ds]] <- norm.list[[ds]] - apply(norm.list[[ds]],1,FUN="min")
      }
    }
  }
  
  ## Min-max the normalized matrix if requested
  if  (norm.minmax) {
    for (ds in 1:length(kernel.list)) { 
      norm.list[[ds]] <- (norm.list[[ds]] - apply(norm.list[[ds]],1,FUN="min")) / (apply(norm.list[[ds]],1,FUN="max") - apply(norm.list[[ds]],1,FUN="min"))
    }
  }
  
  ## Initialize matrices
  if (init == "random") {
    Hmat <- list()
    for (ds in 1:length(kernel.list)) { 
      Hmat[[ds]] <- matrix(runif(nrow(kernel.list[[ds]])*k), nrow = k) 
      colnames(Hmat[[ds]]) <- rownames(kernel.list[[ds]])
    }
    Fmat <- list()
    for (ds in 1:length(kernel.list)) { 
      Fmat[[ds]] <- matrix(runif(nrow(kernel.list[[ds]])*k), ncol = k)
      rownames(Fmat[[ds]]) <- rownames(kernel.list[[ds]])
    }
  } else if (init == "clustering") {
    Hmat <- list()
    for (ds in 1:length(kernel.list)) {
      Hmat[[ds]] <- t(e1071::cmeans(cpca.result$cpca[[ds]], center = k, m = m)$membership)
      colnames(Hmat[[ds]]) <- rownames(kernel.list[[ds]])
    }
    Fmat <- list()
    for (ds in 1:length(kernel.list)) {
      linear <- lm.fit(y = t(kernel.list[[ds]]), x = t(Hmat[[ds]]))
      coefs <- coef(linear)
      coefs[coefs < 0] <- 0
      coefs[is.na(coefs)] <- 0
      Fmat[[ds]] <- t(coefs)
    }
  } else if (init == "pca") {
    Hmat <- list()
    for (ds in 1:length(kernel.list)) {
      Hmat[[ds]] <- t(cpca.result$cpca[[ds]][,1:k] - min(cpca.result$cpca[[ds]][,1:k]))
      colnames(Hmat[[ds]]) <- rownames(kernel.list[[ds]])
    }
    Fmat <- list()
    for (ds in 1:length(kernel.list)) {
      linear <- lm.fit(y = t(kernel.list[[ds]]), x = t(Hmat[[ds]]))
      coefs <- coef(linear)
      coefs[coefs < 0] <- 0
      coefs[is.na(coefs)] <- 0
      Fmat[[ds]] <- t(coefs)
    }
  }
  
  ## Calculate remaining matrices
  Vmat <- list()
  DAmat <- list()
  Wmat <- list()
  for (ds in 1:length(kernel.list)) { 
    Vmat[[ds]] <- colSums(snn.list[[ds]])
    DAmat[[ds]] <- diag(Vmat[[ds]])
    linear <- lm.fit(y = t(norm.list[[ds]]), x = t(Hmat[[ds]]))
    coefs <- coef(linear)
    coefs[ coefs < 0 ] <- 0
    coefs[ is.na(coefs) ] <- 0
    Wmat[[ds]] <- t(coefs)
  }
  
  ## Combine into one list
  data.list <- list()
  for (ds in 1:length(kernel.list)) {
    data.list[[ds]] <- list(ds = ds, kernel = kernel.list[[ds]], norm = norm.list[[ds]], rare = rare.list[[ds]], Hmat = Hmat[[ds]], Fmat = Fmat[[ds]], snn = snn.list[[ds]], Vmat = Vmat[[ds]], DAmat = DAmat[[ds]])
  }
  n_ds <- length(data.list)
  list_names <- names(kernel.list)
  
  ## Remove surplus objects from memory
  kernel.list <- NULL
  norm.list <- NULL
  snn.list <- NULL
  DAmat <- NULL
  Hmat <- NULL
  Fmat <- NULL
  Vmat <- NULL
  rare.list <- NULL
  
  ## Share objects
  if (share.objects) {
    data.list <- share(data.list, copyOnWrite = FALSE)
  }
  times <- c()
  
  ## Save all objects
  if (save_all) {
    Hmat.list <- list()
    Fmat.list <- list()
    Wmat.list <- list()
    Hmat.list[[1]] <- t(do.call("cbind", lapply(data.list, FUN = function(x) { x$Hmat } )))
    Fmat.list[[1]] <- do.call("rbind", lapply(data.list, FUN = function(x) { x$Fmat } ))
    Wmat.list[[1]] <- Wmat
  }
  
  ## Solve
  for (iter in 1:iter.max) {
    # Solve new H matrices
    start <- Sys.time()
    iter.result <- future_lapply(data.list, future.seed = TRUE, FUN = function(x) {
      ## Define indices of other samples
      ds <- x$ds
      js <- seq(1, n_ds,1)
      js <- js[-ds]
      
      ## Update H matrix
      # Numerators
      numerator1 <- t(t(JOINTLY:::matDiMult(t(x$Fmat), x$kernel, ncpu)) * (x$rare * alpha))
      numerator2 <- 2 * mu * x$Hmat
      numerator3 <- lambda * JOINTLY:::matDiMult(x$Hmat, x$snn, ncpu)
      numerator4 <- matrix(nrow = k, ncol = ncol(x$Hmat), 0)
      for (js.idx in js) {
        numerator4 <- numerator4 + ((beta * JOINTLY:::matDiMult(t(Wmat[[js.idx]]), x$norm, ncpu)))  + ((beta * JOINTLY:::matTriMult(t(Wmat[[ds]]), Wmat[[ds]], x$Hmat, ncpu)))
      }
      
      # Denominators
      denom1 = t(t(JOINTLY:::matQuadMult(t(x$Fmat), x$kernel, x$Fmat, x$Hmat, ncpu)) * (x$rare * alpha))
      denom1 = denom1 + 2 * mu * JOINTLY:::matTriMult(x$Hmat, t(x$Hmat), x$Hmat, ncpu)
      denom1 = denom1 + lambda * JOINTLY:::matDiMult(x$Hmat, x$DAmat, ncpu)
      denom2 = matrix(nrow = k, ncol = ncol(x$Hmat), 0)
      for (js.idx in js) {
        denom2 = denom2 + (beta * JOINTLY:::matTriMult(t(Wmat[[js.idx]]), Wmat[[js.idx]], x$Hmat, ncpu))
        denom2 = denom2 + (2 * beta * JOINTLY:::matTriMult(t(Wmat[[ds]]), Wmat[[js.idx]], x$Hmat, ncpu))
        denom2 = denom2 + (beta * JOINTLY:::matDiMult(t(Wmat[[ds]]),x$norm, ncpu))
      }
      
      # Final estimate of Hmat
      H_new <- x$Hmat * ((numerator1 + numerator2 + numerator3 + numerator4) / (denom1 + denom2))
      
      # Update F matrix
      F_new <- x$Fmat * (JOINTLY:::matDiMult(x$kernel, t(H_new), ncpu) / JOINTLY:::matQuadMult(x$kernel, x$Fmat, H_new, t(H_new), ncpu))
      
      # Update W matrix
      linear <- lm.fit(y = t(x$norm), x = t(H_new))
      coefs <- coef(linear)
      coefs[ coefs < 0 ] <- 0
      coefs[ is.na(coefs) ] <- 0
      W_new <- t(coefs)
      
      # Return
      return(list(ds = ds, Hmat = H_new, Fmat = F_new, Wmat = W_new))
    })
    
    # Insert new matrices
    for (ds in 1:length(data.list)) {
      data.list[[iter.result[[ds]]$ds]]$Hmat <- iter.result[[ds]]$Hmat
      data.list[[iter.result[[ds]]$ds]]$Fmat <- iter.result[[ds]]$Fmat
      Wmat[[iter.result[[ds]]$ds]] <- iter.result[[ds]]$Wmat
    }
    
    ## Save all objects
    if (save_all) {
      Hmat.list[[iter+1]] <- t(do.call("cbind", lapply(data.list, FUN = function(x) { x$Hmat } )))
      Fmat.list[[iter+1]] <- do.call("rbind", lapply(data.list, FUN = function(x) { x$Fmat } ))
      Wmat.list[[1]] <- Wmat
    }
    
    # Calculate approximate remaining time
    end <- Sys.time()
    times <- c(times, difftime(end, start, units = "secs"))
    avg.time <- mean(times)
    total.time <- (iter.max - iter) * avg.time
    hours = round(total.time %/% (60 * 60))
    minutes = round((total.time - hours*60*60) %/% 60)
    seconds = round(total.time - hours*60*60 - minutes*60)
    hours_str = ifelse(hours == 0, "", paste0(hours, "H "))
    minutes_str = ifelse((minutes == 0 & hours == 0), "", paste0(minutes, "M "))
    seconds_str = paste0(seconds, "S")
    final_str = paste0(hours_str, minutes_str, seconds_str)
    
    # Calculate percent done
    percent <- iter / iter.max * 100
    
    # Print progress bar
    if (progressbar) {
      if (iter < iter.max) {
        cat(paste(sprintf('\r[%-50s] %d%%', paste(rep('=', percent / 2), collapse = ''), floor(percent)), " - Time remaining: ", final_str, "         ", sep=""))
      } else if (iter == iter.max) { 
        cat(paste(sprintf('\r[%-50s] %d%%', paste(rep('=', percent / 2), collapse = ''), floor(percent)), "                                           ", sep=""))
        cat('\n') }
    }
  }
  
  # Extract final matrices
  Hmat.raw <- lapply(data.list, FUN = function(x) { x$Hmat })
  Fmat.raw <- lapply(data.list, FUN = function(x) { x$Fmat })
  names(Hmat.raw) <- list_names
  names(Fmat.raw) <- list_names
  names(Wmat) <- list_names

  # Scale  
  res <- t(do.call("cbind", Hmat.raw))
  res <- scale(res)
  res <- t(scale(t(res)))
  colnames(res) <- paste("jointly_", 1:ncol(res), sep="")
  
  # Create output
  if (save_all) {
    output <- list(Hmat.scaled = res, Hmat = Hmat.raw, Fmat = Fmat.raw, Wmat = Wmat, allH = Hmat.list, allF = Fmat.list, allW = Wmat.list)
  } else {
    output <- list(Hmat.scaled = res, Hmat = Hmat.raw, Fmat = Fmat.raw, Wmat = Wmat)
  }
  
  # Return
  return(output)
}
