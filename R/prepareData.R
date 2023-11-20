#' @title prepareData
#'
#' @description Function to compute the alpha decay kernels, rareity scores and SNN graphs
#'
#' @param dataset.list A list (per batch) of CPCAs. See \link{cpca}
#' @param k.decay The number of neighbors to use for the decay function [default = 5]
#' @param k.rare The index of the nearest neighbor to use for rare scoring [default = 20]
#' @param k.snn The number of neighbors to use for the SNN calculations [default = 30]
#' @param alpha The power of the decay function [default = 5] 
#' @param threshold Threshold the minimum similarity between cells [default = 1e-4]
#' @param prune The minimum number of normalized neighbors between cells to retain the eddge [default = 1/15]
#' @param rare A boolean (TRUE or FALSE) indicating if rare cell index should be calculated from the distance to the Kth neighbour [default = FALSE]
#' @param verbose Boolean (TRUE or FALSE) indicating to print messages [default = TRUE]
#'
#' @return List (per-batch) of alpha decay kernels.
#' @export
#' @import Matrix
#' @import Seurat

prepareData = function(dataset.list, k.decay = 5, k.rare = 20, k.snn = 30, alpha = 5, threshold = 1e-4, prune = 1/15, rare = FALSE, verbose = TRUE) {
  # TODO: Check that alpha is OK
  # Setup to capture results
  if (verbose) { 
    if (rare) { 
      message("Computing decay kernels, SNN graphs and rareity scores.")
    } else {
      message("Computing decay kernels and SNN graphs.")
    }
  }
  kernel.list <- list()
  rare.list <- list()
  snn.list <- list()
  
  # Loop across datasets and compute the kernel
  for (ds in 1:length(dataset.list)) {
    # Define variables
    data <- dataset.list[[ds]]
    
    # Create pair-wise distance matrix
    distances <- cdist(data)
    rownames(distances) <- colnames(distances) <- rownames(data)
    
    # Get distances for k.decay and k.rare
    selected.distances <- t(apply(distances,1, FUN = function(x) { sort(x, method="quick")[c((k.decay+1), (k.rare+1))] }))
    
    # Calculate rarity score
    if (rare) {
      rare.list[[ds]] <- 1-(1/selected.distances[,2])
    } else {
      rare.list[[ds]] <- rep(1, nrow(data))
    }
    names(rare.list)[ds] <- names(dataset.list)[ds]
    
    # Transform distances by bandwidth
    distances.transformed <- t(t(distances) / selected.distances[,1])
    
    # Calculate distances using decay
    K <- exp(-1 * distances.transformed^alpha)
    
    # Filter distances below threshold
    K[ K < threshold ] <- 0
    K <- K + t(K)
    
    # Insert into list
    kernel.list[[ds]] <- as.matrix(K)
    names(kernel.list)[ds] <- names(dataset.list)[ds]
    
    # Compute SNN
    snn.list[[ds]] <- Seurat::FindNeighbors(data, verbose = FALSE, compute.SNN = TRUE, k.param = k.snn, prune.SNN = prune)$snn
    if (class(snn.list[[ds]])[1] == "Graph") { snn.list[[ds]] <- Matrix::sparseMatrix(i = snn.list[[ds]]@i, p = snn.list[[ds]]@p, x = snn.list[[ds]]@x, dims = snn.list[[ds]]@Dim, dimnames = snn.list[[ds]]@Dimnames, index1 = FALSE) }
    names(snn.list)[ds] <- names(dataset.list)[ds]
  }
  
  # Return
  if (verbose) { message("\n")}
  return(list(kernels = kernel.list, rareity = rare.list, snn = snn.list))
}
