#' @title computeSNN
#'
#' @description Function to compute and prune an SNN graph
#'
#' @param dataset.list A list (per batch) of CPCAs. See \link{cpca}
#' @param k The number of neighbors to use for building the SNN graph [default = "rice" (= (2 * n^(1/3)))]
#' @param prune The minimum number of normalized neighbors between cells to retain the eddge [default = 1/15]
#' @param add The number of neighbors to add in addition to k in neighbourhood search [default = 20]
#' @param rare A boolean (TRUE or FALSE) indicating if rare cell index should be calculated from the distance to the Kth neighbour [default = TRUE]
#'
#' @return A list containing a list (per batch) of SNN graphs and a list (per-batch) of rare cell indices 
#' @export
#' @import Matrix
#' @import RANN

computeSNN <- function(dataset.list, k = "rice", prune = 1/15, add = 20, rare = TRUE) {
  # Setup to capture results
  snn.list <- list()
  rare.list <- list()
  
  # Loop across datasets and compute the kernel
  for (ds in 1:length(dataset.list)) {
    # Calculate k using Rice rule
    if (k == "rice") {
      k.use <- ceiling(2 * (nrow(dataset.list[[ds]])^(1/3)))
    } else {
      k.use <- k
    }
    
    # Get nearest neighbors 
    nn <- RANN::nn2(dataset.list[[ds]], k = k.use+add+1)
    if (rare) {
      rare.list[[ds]] <- 1-(1/nn$nn.dists[,k.use+1])
    } else {
      rare.list[[ds]] <- rep(1, nrow(nn$nn.idx))
    }
    nn <- nn$nn.idx
    nn <- nn[,-1]
    
    # Calculate and setup the SNN
    snn <- sapply(X = 1:nrow(nn), FUN = function(x) { vec <- apply(nn[-x,], y = x, 1, FUN = function(x, y) { sum(x %in% nn[y,]) }); append(vec, 0, after = x-1) })
    colnames(snn) <- rownames(dataset.list[[ds]])
    rownames(snn) <- rownames(dataset.list[[ds]])
    snn <- snn / ((k.use+add) + ((k.use+add) - snn))
    
    # Prune the graph
    snn[ snn < prune] <- 0
    snn <- Matrix::Matrix(snn, sparse = TRUE)
    
    # Save the snn
    snn.list[[ds]] <- snn
  }
  
  # Return
  return(list(snn = snn.list, rare = rare.list))
}