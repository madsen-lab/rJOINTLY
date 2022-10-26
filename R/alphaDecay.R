#' @title alphaDecay
#'
#' @description Function to compute the alpha decay kernels
#'
#' @param data A matrix of CPCAs. See \link{cpca}
#' @param k The number of neighbors to use for the decay function [default = 5]
#' @param alpha The power of the decay function [default = 2] 
#' @param theta A number modifying the bandwidth of the decay function [default = 1e-4]
#'
#' @return A matrix containing the alpha decay kernel
#' @export
#' @import Matrix
#' @import RANN

alphaDecay = function(data, k = 5, alpha = 1, theta = 1e-4) {
  # TODO: Check that alpha is OK
  # Define variables
  k_knn <- k * 20
  bth <- (-log(theta))^(1/alpha)
  N <- nrow(data)
    
  # First round kNN
  nn <- RANN::nn2(data, data, k = k_knn)
  nn.dist <- nn$nn.dists
  nn.idx <- nn$nn.idx
  epsilon <- nn.dist[,k+1]
  below_threshold <- nn.dist[,ncol(nn.dist)] >= bth * epsilon
  idx_thresh <- which(below_threshold)
  if (length(idx_thresh) > 0) {
    K <- exp(-(nn.dist[idx_thresh,] / matrix(rep(epsilon[idx_thresh], ncol(nn.dist)), nrow = length(idx_thresh)))^alpha)
    K[ K <= theta ] <- 0
    K <- as.numeric(K)
    i <- matrix(rep(idx_thresh, ncol(nn.idx)), nrow = length(idx_thresh))
    i <- as.numeric(i)
    idx_tmp <- nn.idx[idx_thresh,]
    j <- as.numeric(idx_tmp)
  } else {
    K <- c()
    i <- c()
    j <- c()
  }
   
  # Expand K till 90% of the data is covered
  while (length(idx_thresh) <= 0.9 * N) {
    k_knn <- min(c(20*k_knn,N))
    data_above <- data[ !(seq(1, nrow(data), 1) %in% idx_thresh),,drop=FALSE]
    epsilon_above <- epsilon[ !(seq(1, length(epsilon), 1) %in% idx_thresh)]
    nn_above <- RANN::nn2(data, data_above, k = k_knn)
    nn_above.dist <- nn_above$nn.dists
    nn_above.idx <- nn_above$nn.idx
    below_threshold_above <- nn_above.dist[,ncol(nn_above.dist)] >= bth * epsilon_above
    idx_thresh_above <- which(below_threshold_above)
    if (length(idx_thresh_above) > 0) {
      K_above <- exp(-(nn_above.dist[idx_thresh_above,] / matrix(rep(epsilon_above[idx_thresh_above], ncol(nn_above.dist)), nrow = length(idx_thresh_above)))^alpha)
      K_above[ K_above <= theta ] <- 0
      K_above <- as.numeric(K_above)
      K <- c(K, K_above)
      idx_above <- which(!below_threshold)
      i_above <- matrix(rep(idx_above[ idx_thresh_above], ncol(nn_above.idx)), nrow = length(idx_thresh_above))
      i_above <- as.numeric(i_above)
      i <- c(i, i_above)
      idx_tmp <- nn_above.idx[idx_thresh_above,]
      j_above <- as.numeric(idx_tmp)
      j <- c(j, j_above)
      below_threshold[ idx_above[idx_thresh_above]] <- TRUE
      idx_thresh <- which(below_threshold)
    }
  }
      
  # Search by radius for the remaining points
  if (length(idx_thresh) < N) {
    data_above <- data[ !(seq(1, nrow(data), 1) %in% idx_thresh),,drop=FALSE]
    epsilon_above <- epsilon[ !(seq(1, length(epsilon), 1) %in% idx_thresh)]
    nn_above <- RANN::nn2(data, data_above, searchtype =  "radius", radius = max(epsilon_above) * bth)
    nn_above.dist <- nn_above$nn.dists
    nn_above.idx <- nn_above$nn.idx
    idx_above <- which(!below_threshold) 
    for (m in 1:nrow(nn_above.idx)) {
      i <- c(i, idx_above[m] * rep(1, length(nn_above.idx[m,])))
      K_above <- exp(-(nn_above.dist[m,] / epsilon_above[m])^alpha)
      K_above[ K_above <= theta ] <- 0
      K <- c(K, K_above)
      j <- c(j, nn_above.idx[m,])
    }
  }
      
  # Create the matrix
  kernel <- Matrix::sparseMatrix(i = i, j = j, x = K)
  kernel <- kernel + Matrix::t(kernel)
    
  # Return
  return(kernel)
}
