#' @title JOINTsolve
#'
#' @description Function to solve the NMF problem in JOINTLY
#' 
#' @param kernel.list A list (per-batch) of kernels. See \link{alphaDecay}
#' @param snn.list A list (per-batch) of SNN graphs. See \link{computeSNN}
#' @param rare.list A list (per-batch) of rare cell indices. See \link{computeSNN}
#' @param cpca.result The results from CPCA analysis. See \link{cpca}
#' @param init The name of the method to use for initialization ("random" or "clustering"). [default = "clustering"]
#' @param k The number of factors to calculate [default = 15]
#' @param m The fuzzy parameter of for cluster initialization [default = 4.5]
#' @param iter.max The maximum number of iterations to perform [default = 100]
#' @param alpha Alpha parameter of the loss function [default = 10]
#' @param mu Mu parameter of the loss function [default = 10]
#' @param lambda Lambda parameter of the loss function [default = 5]
#' @param beta Beta parameter of the loss function [default = 10]
#' @param progressbar A logical indicating whether or not to print a progress bar [default = TRUE]
#' @param ncpu The number of cpus to use for matrix multiplication [default = 1]
#' @param bpparam The name of the type of BPPARAM object to use for sequential or parallel processing [defult = SerialParam()]
#'
#' @return A list containing normalized H matrices to use for clustering, as well a list (per-batch) of H matrices, a list (per-batch) of F matrices and a list (per-batch) of W matrices.
#' @import BiocParallel
#' @import e1071
#' @import SharedObject
#' @export

JOINTLYsolve <- function(kernel.list, snn.list, rare.list, cpca.result, init = "clustering", k = 15, m = 4.5, iter.max = 100, alpha = 10, mu = 10, lambda = 5, beta = 10, progressbar = TRUE, ncpu = 1, bpparam = SerialParam()) {
  ## Convert to dense matrices
  norm.list <- list()
  for (ds in 1:length(kernel.list)) { 
    kernel.list[[ds]] <- as.matrix(kernel.list[[ds]]) 
    snn.list[[ds]] <- as.matrix(snn.list[[ds]]) 
    norm.list[[ds]] <- as.matrix(cpca.result$normalized[[ds]]) 
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
  
  ## Share objects
  kernel.shared <- share(kernel.list)
  norm.shared <- share(norm.list)
  snn.shared <- share(snn.list)
  DA.shared <- share(DAmat)
  times <- c()
  
  ## Solve
  for (iter in 1:iter.max) {
    # Solve new H matrices
    start <- Sys.time()
    iter.result <- bpmapply(function(ds) {
      ## Define indices of other samples
      js <- seq(1, length(kernel.shared),1)
      js <- js[-ds]
      
      ## Update H matrix
      # Numerators
      numerator1 <- t(t(matDiMult(t(Fmat[[ds]]), kernel.shared[[ds]], ncpu)) * (rare.list[[ds]] * alpha))
      numerator2 <- 2 * mu * Hmat[[ds]]
      numerator3 <- lambda * matDiMult(Hmat[[ds]], snn.shared[[ds]], ncpu)
      numerator4 <- matrix(nrow = k, ncol = ncol(Hmat[[ds]]), 0)
      for (js.idx in js) {
        numerator4 <- numerator4 + ((beta * matDiMult(t(Wmat[[js.idx]]), norm.shared[[ds]], ncpu)))  + ((beta * matTriMult(t(Wmat[[ds]]), Wmat[[ds]], Hmat[[ds]], ncpu)))
      }
      
      # Denominators
      denom1 = t(t(matQuadMult(t(Fmat[[ds]]), kernel.shared[[ds]], Fmat[[ds]], Hmat[[ds]], ncpu)) * (rare.list[[ds]] * alpha))
      denom1 = denom1 + 2 * mu * matTriMult(Hmat[[ds]], t(Hmat[[ds]]), Hmat[[ds]], ncpu)
      denom1 = denom1 + lambda * matDiMult(Hmat[[ds]], DA.shared[[ds]], ncpu)
      denom2 = matrix(nrow = k, ncol = ncol(Hmat[[ds]]), 0)
      for (js.idx in js) {
        denom2 = denom2 + (beta * matTriMult(t(Wmat[[js.idx]]), Wmat[[js.idx]], Hmat[[ds]], ncpu))
        denom2 = denom2 + (2 * beta * matTriMult(t(Wmat[[ds]]), Wmat[[js.idx]], Hmat[[ds]], ncpu))
        denom2 = denom2 + (beta * matDiMult(t(Wmat[[ds]]),norm.shared[[ds]], ncpu))
      }
      
      # Final estimate of Hmat
      H_new <- Hmat[[ds]] * ((numerator1 + numerator2 + numerator3 + numerator4) / (denom1 + denom2))
      
      # Update F matrix
      F_new <- Fmat[[ds]] * (matDiMult(kernel.shared[[ds]], t(H_new), ncpu) / matQuadMult(kernel.shared[[ds]], Fmat[[ds]], H_new, t(H_new), ncpu))
      
      # Update W matrix
      linear <- lm.fit(y = t(norm.shared[[ds]]), x = t(H_new))
      coefs <- coef(linear)
      coefs[ coefs < 0 ] <- 0
      coefs[ is.na(coefs) ] <- 0
      W_new <- t(coefs)
      
      # Return
      return(list(Hnew = H_new, Fnew = F_new, Wnew = W_new))
    },1:length(kernel.list), SIMPLIFY = FALSE, BPPARAM = bpparam)
    
    # Insert new matrices
    for (ds in 1:length(kernel.list)) {
      Hmat[[ds]] <- iter.result[[ds]]$Hnew
      Fmat[[ds]] <- iter.result[[ds]]$Fnew
      Wmat[[ds]] <- iter.result[[ds]]$Wnew
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
  
  # Set names
  names(Hmat) <- names(kernel.list)
  names(Fmat) <- names(kernel.list)
  names(Wmat) <- names(kernel.list)
  
  # Scale and combine Hmatrices
  for (i in 1:length(Hmat)) {
    Hmat[[i]] <- t(scale(t(Hmat[[i]])))
    rownames(Hmat[[i]]) <- paste("JOINTLY_", seq(1, k, 1), sep = "")
  }
  res <- t(do.call("cbind", Hmat))
  
  # Return
  return(list(Hmat.scaled = res, Hmat = Hmat, Fmat = Fmat, Wmat = Wmat))
}
