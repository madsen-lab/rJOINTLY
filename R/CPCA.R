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
#' @param do.scale Boolean indicating whether or not to scale expression values prior to decomposition [default = TRUE]
#' @param feat.type The method (either 'inclusive' or 'exclusive') for selecting variable feature [default = "inclusive"]
#' @param oversampling The number of PCs to calculate in excess the selected PCs [default = 10]
#' @param ncpu The number of cpus to use for matrix multiplications [default = 1]
#' @param iter.max The maximum number of iterations to run svd for [default = 100].
#' @param tol Tolerance (minimum increase in fit) for early stopping svd decomposition [default = 1e-5]
#' @param eps Epsilon (minimum values of distances) for early stopping svd decomposition [default = .Machine$double.eps ^ (4 / 5)]
#' @param verbose Boolean (TRUE or FALSE) indicating to print messages [default = TRUE]
#' @param norm.method The method (either 'logCPM' or 'shiftedLog') for normalizing counts [default = 'logCPM']
#' @param sf.method The method (either 'simple', 'scran', 'scranCluster') for normalizing counts [default = 'simple']
#' @param bpparam *Param to use for parallel processing [default = SerialParam()]
#'
#' @return List containing a list (per-batch) of CPCAs and a list (per-batch) of normalized counts.
#' @export
#' @import Seurat
#' @import SingleCellExperiment
#' @import scry
#' @import irlba
#' @import BiocParallel
#' @import scran
#' @import transformGamPoi

cpca = function (dataset.list, weight_by_var = TRUE, pca.type = "cpca", nfeat = 1000, selection.method = "deviance", 
                 threshold = 0.8, kc = 20, ki = 20, do.scale = TRUE, do.center = TRUE, feat.type = "inclusive", oversampling = 10, ncpu = 1, 
                 iter.max = 100, tol = 1e-05, eps = .Machine$double.eps^(4/5), 
                 verbose = TRUE, norm.method = "logCPM", sf.method = "simple", bpparam = SerialParam()) 
{
  if (verbose) { message("Decomposing samples using (consensus) PCA.")}
  if (verbose) { message("\tPerforming feature selection.")}
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
  if (verbose) { message("\tNormalizing and scaling counts.")}
  norm.list <- list()
  for (ds in 1:length(dataset.list)) {
    x <- dataset.list[[ds]]
    
    # Size factors
    if (sf.method == "simple") {
      sf <- Matrix::colSums(x)/10000
    } else if (sf.method == "scranCluster") {
      cl <- scran::quickCluster(x)
      sf <- scran::calculateSumFactors(x, clusters = cl)
      zero.ids <- is.nan(sf) | sf <= 0
      sf[zero.ids] <- NA
      if (any(zero.ids)) {
        sf <- sf/exp(mean(log(sf), na.rm = TRUE))
        sf[zero.ids] <- 0.001
      } else {
      sf <- sf/exp(mean(log(sf)))
      }
    } else if (sf.method == "scran") {
      sf <- scran::calculateSumFactors(x)
      zero.ids <- is.nan(sf) | sf <= 0
      sf[zero.ids] <- NA
      if (any(zero.ids)) {
        sf <- sf/exp(mean(log(sf), na.rm = TRUE))
        sf[zero.ids] <- 0.001
      } else {
      sf <- sf/exp(mean(log(sf)))
      }
    }
    
    # Transformation
    if (norm.method == "logCPM") {
      x <- Matrix::t(Matrix::t(x) / sf)
      x@x <- log1p(x@x)
      x <- x[rownames(x) %in% sel.features, ]
      x <- x[match(sel.features, rownames(x)), ]
      norm.list[[ds]] <- x
    } else if (norm.method == "shiftedLog") {
      x <- transformGamPoi::shifted_log_transform(x, size_factors = sf)
      x <- x[rownames(x) %in% sel.features, ]
      x <- x[match(sel.features, rownames(x)), ]
      norm.list[[ds]] <- x
    }
  }
  scale.list <- list()
  for (ds in 1:length(norm.list)) {
    x <- norm.list[[ds]]
    if (do.scale | do.center) {
      x <- scale(Matrix::t(x), scale = do.scale, center = do.center)
    }
    scale.list[[ds]] <- x
  }
  
  if (pca.type == "cpca") {
    if (verbose) { message("\tCalculating variance-covariance matrix.")}
    V.list <- bplapply(scale.list, BPPARAM = bpparam, FUN = function(x) {
      return(matDiMult(Matrix::t(x), x, n_cores = ncpu)/(nrow(x) - 1))
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
    if (verbose) { message("\tRunning SVD.")}
    for (iter in 1:iter.max) {
      Q <- qr.Q(qr(matDiMult(V, Q, n_cores = ncpu)))
      B <- matDiMult(t(V), Q, n_cores = ncpu)
      Q <- qr.Q(qr(B))
      d_new <- svd(B, nu = 0, nv = 0)$d[1:k]
      idx <- d_new > eps
      if (all(!idx)) 
        break
      if (max(abs((d_new[idx] - d[idx])/d[idx])) < tol) 
        break
      d <- d_new
    }
    Q <- qr.Q(qr(matDiMult(V, Q, n_cores = ncpu)))
    B <- matDiMult(t(Q), V, n_cores = ncpu)
    sc <- svd(B)
    sc$u <- matDiMult(Q, sc$u, n_cores = ncpu)
    sc$u <- sc$u[, 1:k]
    sc$d <- sc$d[1:k]
    sc$v <- sc$v[, 1:k]
    sc$mprod <- 2 * iter + 1
    if (verbose) { message("\tCalculating explained variance.")}
    VScale <- list()
    for (ds in 1:length(scale.list)) {
      VScale[[ds]] <- list(ds = ds, V = V.list[[ds]], S = scale.list[[ds]])
    }
    var.explain <- bplapply(VScale, BPPARAM = bpparam, FUN = function(x) {
      k <- kc
      n <- min(ncol(x$V), k + oversampling)
      Q <- matrix(rnorm(ncol(x$V) * n), ncol(x$V))
      d <- rep(0, k)
      for (iter in 1:iter.max) {
        Q <- qr.Q(qr(matDiMult(x$V, Q, n_cores = ncpu)))
        B <- matDiMult(t(x$V), Q, n_cores = ncpu)
        Q <- qr.Q(qr(B))
        d_new <- svd(B, nu = 0, nv = 0)$d[1:k]
        idx <- d_new > eps
        if (all(!idx)) 
          break
        if (max(abs((d_new[idx] - d[idx])/d[idx])) < tol) 
          break
        d <- d_new
      }
      Q <- qr.Q(qr(matDiMult(x$V, Q, n_cores = ncpu)))
      B <- matDiMult(t(Q), x$V, n_cores = ncpu)
      s <- svd(B)
      s$u <- matDiMult(Q, s$u, n_cores = ncpu)
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
    
    
    if (length(ds.individual) > 0) {
      if (verbose) { message(paste("\tProcessing ", length(ds.individual), " batches with excess variance.", sep=""))}
      SV.list <- list()
      counter <- 1
      for (ds in ds.individual) {
        SV.list[[counter]] <- list(V = V.list[[ds]], S = scale.list[[ds]], ds = ds)
        counter <- counter + 1
      }
      S.list <- bplapply(SV.list, BPPARAM = bpparam, FUN = function(x) {
        V <- x$V - matTriMult(sc$u, t(sc$u), x$V, n_cores = ncpu)
        k <- ki
        n <- min(ncol(V), k + oversampling)
        Q <- matrix(rnorm(ncol(V) * n), ncol(V))
        d <- rep(0, k)
        for (iter in 1:iter.max) {
          Q <- qr.Q(qr(matDiMult(x$V, Q, n_cores = ncpu)))
          B <- matDiMult(t(x$V), Q, n_cores = ncpu)
          Q <- qr.Q(qr(B))
          d_new <- svd(B, nu = 0, nv = 0)$d[1:k]
          idx <- d_new > eps
          if (all(!idx)) 
            break
          if (max(abs((d_new[idx] - d[idx])/d[idx])) < tol) 
            break
          d <- d_new
        }
        Q <- qr.Q(qr(matDiMult(x$V, Q, n_cores = ncpu)))
        B <- matDiMult(t(Q), x$V, n_cores = ncpu)
        s <- svd(B)
        s$u <- matDiMult(Q, s$u, n_cores = ncpu)
        s$u <- s$u[, 1:k]
        s$d <- s$d[1:k]
        s$v <- s$v[, 1:k]
        s$mprod <- 2 * iter + 1
        var.exp <- c()
        total <- t(matDiMult(t(sc$u), t(x$S), n_cores = ncpu))
        for (k.sel in 1:ki) {
          var.exp <- c(var.exp, sum(apply(cbind(total, t(matDiMult(t(s$u[, 1:k.sel, drop = FALSE]), t(x$S), n_cores = ncpu))), 2, FUN = "var"))/sum(apply(x$S, 2, FUN = "var")))
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
      C.list[[ds]] <- t(matDiMult(t(common), t(scale.list[[ds]]), n_cores = ncpu))
      rownames(C.list[[ds]]) <- colnames(dataset.list[[ds]])
      names(C.list)[ds] <- names(dataset.list)[ds]
    }
  } else {
    if (verbose) { message("\tRunning per-batch SVD.")}
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
  if (verbose) { message("\n")}
  return(list(cpca = C.list, normalized = norm.list, features = sel.features))
}
