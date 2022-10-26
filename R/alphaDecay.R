#' @title alphaDecay
#'
#' @description Function to compute the alpha decay kernels and rareity scores
#'
#' @param dataset.list A list (per batch) of CPCAs. See \link{cpca}
#' @param k.decay The number of neighbors to use for the decay function [default = 5]
#' @param k.rare The index of the nearest neighbor to use for rare scoring [default = "rice" (= (2 * n^(1/3)))]
#' @param alpha The power of the decay function [default = 1] 
#' @param threshold Threshold the minimum similarity between cells [default = 1e-4]
#' @param rare A boolean (TRUE or FALSE) indicating if rare cell index should be calculated from the distance to the Kth neighbour [default = TRUE]
#'
#' @return List (per-batch) of alpha decay kernels.
#' @export
#' @import Matrix

alphaDecay = function(dataset.list, k.decay = 5, k.rare = "rice", alpha = 1, threshold = 1e-4, rare = TRUE) {
  # TODO: Check that alpha is OK
  # Setup to capture results
  kernel.list <- list()
  rare.list <- list()
  
  # Loop across datasets and compute the kernel
  for (ds in 1:length(dataset.list)) {
    # Define variables
    data <- dataset.list[[ds]]
    
    # Create pair-wise distance matrix
    distances <- cdist(data)
    
    # Calculate rarity score
    if (k.rare == "rice") {
      k.use <- ceiling(2 * (nrow(data)^(1/3)))
    } else {
      k.use <- k
    }
    if (rare) {
      rare.list[[ds]] <- 1-(1/apply(distances,1, FUN = function(x) { sort(x, method="quick")[k.use+1]} ))
    } else {
      rare.list[[ds]] <- rep(1, nrow(data))
    }
    
    # Calculate bandwidth from the maximum distance to the k+1 nearest neighbor 
    bandwidth <- apply(distances,1, FUN = function(x) { sort(x, method="quick")[k.decay+1]} )
    
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
  return(list(kernels = kernel.list, rareity = rare.list))
}
