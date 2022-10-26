#' @title alphaDecay
#'
#' @description Function to compute the alpha decay kernels
#'
#' @param dataset.list A list (per batch) of CPCAs. See \link{cpca}
#' @param k The number of neighbors to use for the decay function [default = 5]
#' @param alpha The power of the decay function [default = 2] 
#' @param threshold Threshold the minimum similarity between cells [default = 1e-4]
#'
#' @return List (per-batch) of alpha decay kernels.
#' @export
#' @import Matrix

alphaDecay = function(dataset.list, k = 5, alpha = 1, threshold = 1e-4) {
  # TODO: Check that alpha is OK
  # Setup to capture results
  kernel.list <- list()
  
  # Loop across datasets and compute the kernel
  for (ds in 1:length(dataset.list)) {
    # Define variables
    data <- dataset.list[[ds]]
    
    # Create pair-wise distance matrix
    distances <- cdist(data)
    
    # Calculate bandwidth from the maximum distance to the k+1 nearest neighbor 
    bandwidth <- apply(distances,1, FUN = function(x) { sort(x, method="quick")[k+1]} )
    
    # Transform distances by bandwidth
    distances.transformed <- t(t(distances) / bandwidth)
    
    # Calculate distances using decay
    K <- exp(-1 * distances.transformed^alpha)
    
    # Filter distances below threshold
    K[ K < threshold ] <- 0
    K <- K + t(K)
      
    # Insert into list
    kernel.list[[ds]] <- as.matrix(K)
    names(kernel.list)[ds] <- names(dataset.list)[ds]
  }
  
  # Return
  return(kernel.list)
}
