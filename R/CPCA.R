#' @title CPCA
#'
#' @description Function to perform consensus PCA
#'
#' @param dataset.list A list (per batch) of raw counts
#' @param weight_by_var Boolean (TRUE or FALSE) indicating whether or not to weight embeddings by the variance of each PC [default = TRUE]
#' @param pca.type The method (either 'cpca' or 'irlba') for computing PCs [default = "cpca"]
#' @param nfeat The number of features to select [default = 1000]
#' @param selection.method The method to use selecting highly-variable features [default = "deviance"]
#' @param threshold The minimum amount of variance to be explained by common PCs before adding individual PCs [default = 0.8]
#' @param kc The number of common PCs to calculate [default = 20]
#' @param ki The number of invididual PCs to calcuate [default = 20]
#' @param do.center Boolean indicating whether or not to center expression values prior to decomposition [default = TRUE]
#' @param feat.type The method (either 'inclusive' or 'exclusive') for selecting variable feature [default = "inclusive"]
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
#' @import irlba

cpca = function (dataset.list, weight_by_var = TRUE, pca.type = "cpca", nfeat = 1000, selection.method = "deviance", 
                 threshold = 0.8, kc = 20, ki = 20, do.center = TRUE, feat.type = "inclusive", oversampling = 10, ncpu = 1, 
                 iter.max = 100, tol = 1e-05, eps = .Machine$double.eps^(4/5), 
                 verbose = TRUE) 
{
  if (verbose) { message("Performing feature selection.")}
  if (selection.method == "deviance") {
    features <- as.data.frame(scry::devianceFeatureSelection(dataset.list[[1]]))
    colnames(features)[1] <- "Dataset1"
    for (ds in 2:length(dataset.list)) {
      feats <- scry::devianceFeatureSelection(dataset.list[[ds]])
      features <- merge(as.data.frame(features), as.data.frame(feats), 
                        by = 0)
      colnames(features)[ncol(features)] <- paste("Dataset", 
                                                  ds, sep = "")
      rownames(features) <- features[, 1]
      features <- features[, -1]
    }
    for (col in 1:ncol(features)) {
      features <- features[!is.na(features[, col]), ]
    }
    features$average <- apply(features[, c(1:ncol(features))], 
                              1, FUN = "mean")
    features$count <- 0
    for (col in 1:(ncol(features) - 2)) {
      features[which(rank(-features[, col]) <= nfeat), 
               "count"] <- features[which(rank(-features[, col]) <= 
                                            nfeat), "count"] + 1
    }
    features <- features[order(features$count, features$average, 
                               decreasing = TRUE), ]
    if (feat.type == "inclusive") {
      sel.features <- rownames(features[features$count > 0, 
      ])
    } else {
      sel.features <- rownames(features[1:nfeat,])
    }
  }
  else if (selection.method == "seurat") {
    seu.list <- list()
    for (ds in 1:length(dataset.list)) {
      seu.list[[ds]] <- Seurat::CreateSeuratObject(counts = dataset.list[[ds]])
    }
    sel.features <- Seurat::SelectIntegrationFeatures(seu.list, 
                                                      nfeatures = nfeat, verbose = FALSE)
  }
  else if (selection.method == "none") {
    sel.features <- intersect(rownames(dataset.list[[1]]), 
                              rownames(dataset.list[[2]]))
    if (length(dataset.list) >= 3) {
      for (ds in 3:length(dataset.list)) {
        sel.features <- intersect(sel.features, rownames(dataset.list[[ds]]))
      }
    }
  }
  else {
  }
  if (length(sel.features) == 0) {
  }
  if (verbose) { message("Normalizing and scaling counts.")}
  norm.list <- list()
  for (ds in 1:length(dataset.list)) {
    x <- dataset.list[[ds]]
    sf <- 10000/Matrix::colSums(x)
    x <- Matrix::t(Matrix::t(x) * sf)
    x@x <- log1p(x@x)
    x <- x[rownames(x) %in% sel.features, ]
    x <- x[match(sel.features, rownames(x)), ]
    norm.list[[ds]] <- x
  }
  scale.list <- list()
  for (ds in 1:length(norm.list)) {
    x <- norm.list[[ds]]
    x <- scale(Matrix::t(x), center = do.center)
    scale.list[[ds]] <- x
  }
  
  if (pca.type == "cpca") {
    if (verbose) { message("Calculating variance-covariance matrix.")}
    V.list <- future_lapply(scale.list, future.seed = TRUE, FUN = function(x) {
      return(JOINTLY:::matDiMult(Matrix::t(x), x, n_cores = ncpu)/(nrow(x) - 1))
    })
    V <- matrix(nrow = dim(V.list[[1]])[1], ncol = dim(V.list[[1]])[1], 0)
    sums <- sum(unlist(lapply(scale.list, FUN = "nrow")))
    for (ds in 1:length(scale.list)) {
      V = V + V.list[[ds]] * nrow(scale.list[[ds]])/sums
    }
    k <- kc
    n <- min(nrow(V), k + oversampling)
    Q <- matrix(rnorm(nrow(V) * n), nrow(V))
    d <- rep(0, k)
    if (verbose) { message("Running consensus PCA.")}
    for (iter in 1:iter.max) {
      Q <- qr.Q(qr(JOINTLY:::matDiMult(V, Q, n_cores = ncpu)))
      B <- JOINTLY:::matDiMult(t(V), Q, n_cores = ncpu)
      Q <- qr.Q(qr(B))
      d_new <- svd(B, nu = 0, nv = 0)$d[1:k]
      idx <- d_new > eps
      if (all(!idx)) 
        break
      if (max(abs((d_new[idx] - d[idx])/d[idx])) < tol) 
        break
      d <- d_new
    }
    Q <- qr.Q(qr(JOINTLY:::matDiMult(V, Q, n_cores = ncpu)))
    B <- JOINTLY:::matDiMult(t(Q), V, n_cores = ncpu)
    sc <- svd(B)
    sc$u <- JOINTLY:::matDiMult(Q, sc$u, n_cores = ncpu)
    sc$u <- sc$u[, 1:k]
    sc$d <- sc$d[1:k]
    sc$v <- sc$v[, 1:k]
    sc$mprod <- 2 * iter + 1
    if (verbose) { message("Calculating explained variance.")}
    VScale <- list()
    for (ds in 1:length(scale.list)) {
      VScale[[ds]] <- list(ds = ds, V = V.list[[ds]], S = scale.list[[ds]])
    }
    var.explain <- future_lapply(VScale, future.seed = TRUE, FUN = function(x) {
      k <- kc
      n <- min(ncol(x$V), k + oversampling)
      Q <- matrix(rnorm(ncol(x$V) * n), ncol(x$V))
      d <- rep(0, k)
      for (iter in 1:iter.max) {
        Q <- qr.Q(qr(JOINTLY:::matDiMult(x$V, Q, n_cores = ncpu)))
        B <- JOINTLY:::matDiMult(t(x$V), Q, n_cores = ncpu)
        Q <- qr.Q(qr(B))
        d_new <- svd(B, nu = 0, nv = 0)$d[1:k]
        idx <- d_new > eps
        if (all(!idx)) 
          break
        if (max(abs((d_new[idx] - d[idx])/d[idx])) < tol) 
          break
        d <- d_new
      }
      Q <- qr.Q(qr(JOINTLY:::matDiMult(x$V, Q, n_cores = ncpu)))
      B <- JOINTLY:::matDiMult(t(Q), x$V, n_cores = ncpu)
      s <- svd(B)
      s$u <- JOINTLY:::matDiMult(Q, s$u, n_cores = ncpu)
      s$u <- s$u[, 1:k]
      s$d <- s$d[1:k]
      s$v <- s$v[, 1:k]
      s$mprod <- 2 * iter + 1
      return(c(x$ds, sum(apply(t(t(sc$u) %*% t(x$S)), 
                               2, FUN = "var"))/sum(apply(x$S, 2, FUN = "var")), sum(apply(t(t(s$u) %*% t(x$S)), 
                                                                                           2, FUN = "var"))/sum(apply(x$S, 2, FUN = "var"))))
    })
    var.explain <- as.data.frame(do.call("rbind",var.explain))
    ds.individual <- var.explain[var.explain[, 2]/var.explain[, 3] < threshold, 1]
    if (verbose) { message("Processing samples with excess variance.")}
    
    if (length(ds.individual) > 0) {
      SV.list <- list()
      counter <- 1
      for (ds in ds.individual) {
        SV.list[[counter]] <- list(V = V.list[[ds]], S = scale.list[[ds]], ds = ds)
        counter <- counter + 1
      }
      S.list <- future_lapply(SV.list, future.seed = TRUE, FUN = function(x) {
        V <- x$V - JOINTLY:::matTriMult(sc$u, t(sc$u), x$V, n_cores = ncpu)
        k <- ki
        n <- min(ncol(V), k + oversampling)
        Q <- matrix(rnorm(ncol(V) * n), ncol(V))
        d <- rep(0, k)
        for (iter in 1:iter.max) {
          Q <- qr.Q(qr(JOINTLY:::matDiMult(x$V, Q, n_cores = ncpu)))
          B <- JOINTLY:::matDiMult(t(x$V), Q, n_cores = ncpu)
          Q <- qr.Q(qr(B))
          d_new <- svd(B, nu = 0, nv = 0)$d[1:k]
          idx <- d_new > eps
          if (all(!idx)) 
            break
          if (max(abs((d_new[idx] - d[idx])/d[idx])) < tol) 
            break
          d <- d_new
        }
        Q <- qr.Q(qr(JOINTLY:::matDiMult(x$V, Q, n_cores = ncpu)))
        B <- JOINTLY:::matDiMult(t(Q), x$V, n_cores = ncpu)
        s <- svd(B)
        s$u <- JOINTLY:::matDiMult(Q, s$u, n_cores = ncpu)
        s$u <- s$u[, 1:k]
        s$d <- s$d[1:k]
        s$v <- s$v[, 1:k]
        s$mprod <- 2 * iter + 1
        var.exp <- c()
        total <- t(JOINTLY:::matDiMult(t(sc$u), t(x$S), n_cores = ncpu))
        for (k.sel in 1:ki) {
          var.exp <- c(var.exp, sum(apply(cbind(total, t(JOINTLY:::matDiMult(t(s$u[, 1:k.sel, drop = FALSE]), t(x$S), n_cores = ncpu))), 2, FUN = "var"))/sum(apply(x$S, 2, FUN = "var")))
        }
        k.sel <- which.min(abs(var.exp - var.explain[var.explain[, 1] == x$ds, 3] * threshold))
        s$u <- s$u[, 1:k.sel]
        s$d <- s$d[1:k.sel]
        s$v <- s$v[, 1:k.sel]
        return(s)
      })
    }
 
    common <- sc$u
    if (length(ds.individual) > 0) {
      for (ds in 1:length(S.list)) {
        invidual <- S.list[[ds]]$u
        common <- cbind(common, invidual)
      }
    }
    C.list <- list()
    for (ds in 1:length(scale.list)) {
      C.list[[ds]] <- t(JOINTLY:::matDiMult(t(common), t(scale.list[[ds]]), n_cores = ncpu))
      rownames(C.list[[ds]]) <- colnames(dataset.list[[ds]])
      names(C.list)[ds] <- names(dataset.list)[ds]
    }
  } else {
    C.list <- list()
    scale.data <- do.call("rbind", scale.list)
    svd <- irlba::irlba(scale.data, nv = kc)
    if (weight_by_var) {
      embed <- svd$u %*% diag(svd$d)
    } else {
      embed <- svd$u
    }
    rownames(embed) <- rownames(scale.data)
    for (ds in 1:length(scale.list)) {
      C.list[[ds]] <- embed[ rownames(embed) %in% rownames(scale.list[[ds]]),]
    }
  }
  return(list(cpca = C.list, normalized = norm.list, features = sel.features))
}