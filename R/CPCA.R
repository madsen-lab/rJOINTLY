#' @title CPCA
#'
#' @description Function to perform consensus PCA
#'
#' @param dataset.list A list (per batch) of raw counts
#' @param nfeat The number of features to select [default = 1000]
#' @param selection.method The method to use selecting highly-variable features [default = "deviance"]
#' @param threshold The minimum amount of variance to be explained by common PCs before adding individual PCs [default = 0.8]
#' @param kc The number of common PCs to calculate [default = 20]
#' @param ki The number of invididual PCs to calcuate [default = 20]
#' @param oversampling The number of PCs to calculate in excess the selected PCs [default = 10]
#' @param ncpu The number of cpus to use for matrix multiplications [default = 1]
#' @param iter.max The maximum number of iterations to run svd for [default = 100].
#' @param tol Tolerance (minimum increase in fit) for early stopping svd decomposition [default = 1e-5]
#' @param eps Epsilon (minimum values of distances) for early stopping svd decomposition [default = .Machine$double.eps ^ (4 / 5)]
#' @param verbose Boolean (TRUE or FALSE) indicating to print messages [default = TRUE]
#'
#' @return List containing a list (per-batch) of CPCAs and a list (per-batch) of normalized counts.
#' @export
#' @import Seurat
#' @import SingleCellExperiment
#' @import scry

cpca = function(dataset.list, nfeat = 1000, selection.method = "deviance", threshold = 0.8, kc = 20, ki = 20, oversampling = 10, ncpu = 1, iter.max = 100, tol = 1e-5, eps = .Machine$double.eps ^ (4 / 5), verbose = TRUE) {
  ## TODO: Check parameters
  
  ## Feature selection
  ## TO DO: User defined features
  if (selection.method == "deviance") {
    ### Select features
    features <- as.data.frame(scry::devianceFeatureSelection(dataset.list[[1]]))
    colnames(features)[1] <- "Dataset1"
    for (ds in 2:length(dataset.list)) {
      feats <- scry::devianceFeatureSelection(dataset.list[[ds]])
      features <- merge(as.data.frame(features), as.data.frame(feats), by=0)
      colnames(features)[ncol(features)] <- paste("Dataset", ds, sep="")
      rownames(features) <- features[,1]
      features <- features[,-1]
    }
    for (col in 1:ncol(features)) { features <- features[ !is.na(features[,col]),] }
    features$average <- apply(features[,c(1:ncol(features))],1,FUN="mean")
    features$count <- 0
    for (col in 1:(ncol(features)-2)) { features[ which(rank(-features[,col]) <= nfeat),"count"] <- features[ which(rank(-features[,col]) <= nfeat),"count"] + 1 }
    features <- features[ order(features$count, features$average, decreasing = TRUE),]
    sel.features <- rownames(features[ features$count > 0,])
  } else if (selection.method == "seurat") {
    # Select features using Seurat
    seu.list <- list()
    for (ds in 1:length(dataset.list)) {
      seu.list[[ds]] <- Seurat::CreateSeuratObject(counts = dataset.list[[ds]])
    }
    sel.features <- Seurat::SelectIntegrationFeatures(seu.list, nfeatures = nfeat, verbose = FALSE)
  } else if (selection.method == "none") {
    sel.features <- intersect(rownames(dataset.list[[1]]), rownames(dataset.list[[2]]))
    if (length(dataset.list) >= 3) {
      for (ds in 3:length(dataset.list)) { sel.features <- intersect(sel.features, rownames(dataset.list[[ds]])) } 
    }
  } else {
    ## Exit with message
  }
  
  # Check how many variable features were selected
  if (length(sel.features) == 0) {
    ## Break with error massage
  }
  
  ## TODO: Normalize data and subset
  norm.list <-  list()
  for (ds in 1:length(dataset.list)) {
    x <- dataset.list[[ds]]
    sf <- 10000 / Matrix::colSums(x)
    x <- Matrix::t(Matrix::t(x) * sf)
    x@x <- log1p(x@x)
    x <- x[ rownames(x) %in% sel.features,]
    x <- x[ match(sel.features, rownames(x)),]
    norm.list[[ds]] <- x
  }
  
  ## Subset list to selected features and scale
  scale.list <- list()
  for (ds in 1:length(norm.list)) {
    x <- norm.list[[ds]]
    x <- scale(Matrix::t(x))
    scale.list[[ds]] <- x
  }

  ## Calculate variance-covariance matrices
  V.list <- list()
  for (ds in 1:length(scale.list)) {
    V.list[[ds]] <- matDiMult(Matrix::t(scale.list[[ds]]), scale.list[[ds]], n_cores = ncpu) / (nrow(scale.list[[ds]]) - 1)
  }

  ## Within group variance-covariance matrix
  V <- matrix(nrow = dim(V.list[[1]])[1], ncol = dim(V.list[[1]])[1], 0)
  for (ds in 1:length(scale.list)) {
    V = V + V.list[[ds]]*nrow(scale.list[[ds]])/sum(unlist(lapply(scale.list,FUN="nrow")))
  }

  ## Randomized SVD using QR decomposition on the V matrix
  # Setup
  k <- kc
  n <- min(ncol(V), k + oversampling)
  Q <-  matrix(rnorm(ncol(V) * n), ncol(V))
  d <- rep(0, k)

  # Iterate
  for (iter in 1:iter.max) {
    # QR decomposition
    Q <- qr.Q(qr(V %*% Q))
    B <- crossprod(V, Q)
    Q <- qr.Q(qr(B))
    
    # Distance
    d_new <- svd(B, nu=0, nv=0)$d[1:k]
    
    # Evaluate criteria
    idx <- d_new > eps
    if (all(! idx)) break
    if (max(abs((d_new[idx] - d[idx]) / d[idx])) < tol) break
    d <- d_new
  }

  # Finalize results
  Q <- qr.Q(qr(V %*% Q))
  B <- crossprod(Q, V)
  sc <- svd(B)
  sc$u <- Q %*% sc$u
  sc$u <- sc$u[, 1:k]
  sc$d <- sc$d[1:k]
  sc$v <- sc$v[, 1:k]
  sc$mprod <- 2 * iter + 1

  # Check if datasets are integratable
  var.explain <- as.data.frame(matrix(ncol=kc, nrow=length(scale.list)))
  for (ds in 1:length(scale.list)) { var.explain[ds, ] <- apply(t(t(sc$u) %*% t(scale.list[[ds]])),2,FUN="var")/sum(apply(scale.list[[ds]],2,FUN="var")) }

  # Evalute the need for individual components
  var.explain <- as.data.frame(matrix(ncol=3, nrow=length(scale.list)))
  for (ds in 1:length(scale.list)) {
    k <- kc
    V <- V.list[[ds]]
    n <- min(ncol(V), k + oversampling)
    Q <- matrix(rnorm(ncol(V) * n), ncol(V))
    d <- rep(0, k)
    
    # Iterate
    for (iter in 1:iter.max) {
      # QR decomposition
      Q <- qr.Q(qr(V %*% Q))
      B <- crossprod(V, Q)
      Q <- qr.Q(qr(B))
      
      # Distance
      d_new <- svd(B, nu=0, nv=0)$d[1:k]
      
      # Evaluate criteria
      idx <- d_new > eps
      if (all(! idx)) break
      if (max(abs((d_new[idx] - d[idx]) / d[idx])) < tol) break
      d <- d_new
    }
    
    # Finalize results
    Q <- qr.Q(qr(V %*% Q))
    B <- crossprod(Q, V)
    s <- svd(B)
    s$u <- Q %*% s$u
    s$u <- s$u[, 1:k]
    s$d <- s$d[1:k]
    s$v <- s$v[, 1:k]
    s$mprod <- 2 * iter + 1
    
    # Save variance
    var.explain[ds,1] <- ds
    var.explain[ds,2] <- sum(apply(t(t(sc$u) %*% t(scale.list[[ds]])),2,FUN="var"))/sum(apply(scale.list[[ds]],2,FUN="var"))
    var.explain[ds,3] <- sum(apply(t(t(s$u) %*% t(scale.list[[ds]])),2,FUN="var"))/sum(apply(scale.list[[ds]],2,FUN="var"))
  }
  ds.individual <- var.explain[ var.explain[,2] / var.explain[,3] < threshold, 1]

  ## Calculate residuals
  R.list <- structure(vector(mode = "list", length = length(ds.individual)), names = names(dataset.list)[ds.individual])
  for (ds in ds.individual) {
    R.list[[ds]] <- V.list[[ds]] - matTriMult(sc$u, t(sc$u), V.list[[ds]], n_cores = ncpu)
  }

  ## Randomized SVD using QR decomposition on each R matrix
  S.list <- structure(vector(mode = "list", length = length(scale.list)), names = names(dataset.list))
  for (ds in ds.individual) {
    # Setup
    k <- ki
    V <- R.list[[ds]]
    n <- min(ncol(V), k + oversampling)
    Q <-  matrix(rnorm(ncol(V) * n), ncol(V))
    d <- rep(0, k)
  
    # Iterate
    for (iter in 1:iter.max) {
      # QR decomposition
      Q <- qr.Q(qr(V %*% Q))
      B <- crossprod(V, Q)
      Q <- qr.Q(qr(B))
      
      # Distance
      d_new <- svd(B, nu=0, nv=0)$d[1:k]
      
      # Evaluate criteria
      idx <- d_new > eps
      if (all(! idx)) break
      if (max(abs((d_new[idx] - d[idx]) / d[idx])) < tol) break
      d <- d_new
    }
    
    # Finalize results
    Q <- qr.Q(qr(V %*% Q))
    B <- crossprod(Q, V)
    s <- svd(B)
    s$u <- Q %*% s$u
    s$u <- s$u[, 1:k]
    s$d <- s$d[1:k]
    s$v <- s$v[, 1:k]
    s$mprod <- 2 * iter + 1
    
    # Evaluate number of components
    var.exp <- c()
    for (k.sel in 1:ki) {
      var.exp <- c(var.exp, sum(apply(cbind(t(t(sc$u) %*% t(scale.list[[ds]])),t(t(s$u[,1:k.sel,drop=FALSE]) %*% t(scale.list[[ds]]))),2,FUN="var"))/sum(apply(scale.list[[ds]],2,FUN="var")))
    }
    k.sel <- which.min(abs(var.exp - var.explain[ var.explain[,1] == ds, 3] * threshold))
    s$u <- s$u[, 1:k.sel]
    s$d <- s$d[1:k.sel]
    s$v <- s$v[, 1:k.sel]
    S.list[[ds]] <- s
  }
  
  ## Final PC space
  common <- sc$u
  for (ds in 1:length(dataset.list)) {
    if (!is.null(S.list[[ds]])) { 
      invidual <- S.list[[ds]]$u
      common <- cbind(common, invidual)
    }
  }
  
  ## Final rotations
  C.list <- list()
  for (ds in 1:length(scale.list)) {
    C.list[[ds]] <- t(t(common) %*% t(scale.list[[ds]]))
    rownames(C.list[[ds]]) <- colnames(dataset.list[[ds]])
    names(C.list)[ds] <- names(dataset.list)[ds]
  }
  
  ## Return
  return(list(cpca = C.list, normalized = norm.list, svd = sc))
}
