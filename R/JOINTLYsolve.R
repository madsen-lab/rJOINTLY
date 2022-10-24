#' @title JOINTsolve
#'
#' @description Function to solve NMF problems
#' 
#' @param kernel.list A list (per-batch) of kernels. See \link{alphaDecay}
#' @param snn.list A list (per-batch) of SNN graphs. See \link{computeSNN}
#' @param rare.list A list (per-batch) of rare cell indices. See \link{computeSNN}
#' @param norm.list A list (per-batch) of normalized counts for highly variable genes. See \link{cpca}
#' @param k The number of factors to calculate [default = 20]
#' @param iter.max The maximum number of iterations to perform [default = 100]
#' @param alpha Alpha parameter of the loss function [default = 100]
#' @param mu Mu parameter of the loss function [default = 1]
#' @param lambda Lambda parameter of the loss function [default = 100]
#' @param beta Beta parameter of the loss function [default = 1]
#' @param ncpu The number of cpus to use for matrix multiplication [default = 1]
#'
#' @return A list containing a list (per-batch) of H matrices, a list (per-batch) of F matrices and a list (per-batch) of W matrices.
#' @export

JOINTLYsolve <- function(kernel.list, snn.list, rare.list, norm.list, k = 20, iter.max = 100, alpha = 100, mu = 1, lambda = 100, beta = 1, ncpu = 1) {
  ## Convert to dense matrices
  for (ds in 1:length(kernel.list)) { 
    kernel.list[[ds]] <- as.matrix(kernel.list[[ds]]) 
    snn.list[[ds]] <- as.matrix(snn.list[[ds]]) 
    norm.list[[ds]] <- as.matrix(norm.list[[ds]]) 
  }
  
  ## Initialize matrices
  Hmat <- list()
  for (ds in 1:length(kernel.list)) { Hmat[[ds]] <- matrix(runif(nrow(kernel.list[[ds]])*k), nrow = k) }
  Hmat_new <- Hmat
  Fmat <- list()
  for (ds in 1:length(kernel.list)) { Fmat[[ds]] <- matrix(runif(nrow(kernel.list[[ds]])*k), ncol = k) }
  Fmat_new <- Fmat

  ## Calculate matrices
  Vmat <- list()
  for (ds in 1:length(kernel.list)) { Vmat[[ds]] <- colSums(snn.list[[ds]]) }
  DAmat <- list()
  for (ds in 1:length(kernel.list)) { DAmat[[ds]] <- diag(Vmat[[ds]]) }
  Wmat <- list()
  for (ds in 1:length(kernel.list)) {
    linear <- lm.fit(y = t(norm.list[[ds]]), x = t(Hmat[[ds]]))
    coefs <- coef(linear)
    coefs[ coefs < 0 ] <- 0
    Wmat[[ds]] <- t(coefs)
  }
  Wmat_new <- Wmat
  
  ## Solve
  for (iter in 1:iter.max) {
    # Loop across datasets
    for (ds in 1:length(kernel.list)) {
      ## Define indices of other samples
      js <- seq(1, length(kernel.list),1)
      js <- js[-ds]
      
      ## Update H matrix
      # Numerators
      numerator1 <- t(t(matDiMult(t(Fmat[[ds]]), kernel.list[[ds]], ncpu)) * (rare.list[[ds]] * alpha))
      numerator2 <- 2 * mu * Hmat[[ds]]
      numerator3 <- lambda * matDiMult(Hmat[[ds]], snn.list[[ds]], ncpu)
      numerator4 <- matrix(nrow = k, ncol = ncol(Hmat[[ds]]), 0)
      for (js.idx in js) {
        numerator4 <- numerator4 + ((beta * matDiMult(t(Wmat[[js.idx]]), norm.list[[ds]], ncpu)))  + ((beta * matTriMult(t(Wmat[[ds]]), Wmat[[ds]], Hmat[[ds]], ncpu)))
      }
    
      # Denominators
      denom1 = t(t(matQuadMult(t(Fmat[[ds]]), kernel.list[[ds]], Fmat[[ds]], Hmat[[ds]], ncpu)) * (rare.list[[ds]] * alpha))
      denom1 = denom1 + 2 * mu * matTriMult(Hmat[[ds]], t(Hmat[[ds]]), Hmat[[ds]], ncpu) # Wierd to have Hmat twice here?
      denom1 = denom1 + lambda * matDiMult(Hmat[[ds]], DAmat[[ds]], ncpu)
      denom2 = matrix(nrow = k, ncol = ncol(Hmat[[ds]]), 0)
      for (js.idx in js) {
        denom2 = denom2 + (beta * matTriMult(t(Wmat[[js.idx]]), Wmat[[js.idx]], Hmat[[ds]], ncpu))
        denom2 = denom2 + (2 * beta * matTriMult(t(Wmat[[ds]]), Wmat[[js.idx]], Hmat[[ds]], ncpu))
        denom2 = denom2 + (beta * matDiMult(t(Wmat[[ds]]),norm.list[[ds]], ncpu))
      }
      
      # Final estimate of Hmat
      Hmat_new[[ds]] <- Hmat[[ds]] * ((numerator1 + numerator2 + numerator3 + numerator4) / (denom1 + denom2))
      
      ## Update F matrix
      Fmat_new[[ds]] <- Fmat[[ds]] * (matDiMult(kernel.list[[ds]], t(Hmat_new[[ds]]), ncpu) / matQuadMult(kernel.list[[ds]], Fmat[[ds]], Hmat_new[[ds]], t(Hmat_new[[ds]]), ncpu))
      
      ## Update W matrix
      linear <- lm.fit(y = t(norm.list[[ds]]), x = t(Hmat_new[[ds]]))
      coefs <- coef(linear)
      coefs[ coefs < 0 ] <- 0
      Wmat_new[[ds]] <- t(coefs)
    }
    
    # Update all matrices
    Hmat <- Hmat_new
    Fmat <- Fmat_new
    Wmat <- Wmat_new
    print(iter)
  }
  
  # Return
  return(list(Hmat = Hmat, Fmat = Fmat, Wmat = Wmat))
}