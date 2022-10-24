#' @title Preprocess
#'
#' @description Function to preprocess data from Seurat, SCE or Matrix-like objects
#'
#' @param data A single (or list of) Seurat, SingleCellExperiment, or Matrix-like object(s) containing raw counts.
#' @param batch.var A variable to split Seurat or SingleCellExperiment objects into a list [default = NULL]
#'
#' @return A list (per-batch) of raw counts
#' @export
#' @import Seurat
#' @import SingleCellExperiment

preprocess = function(data, batch.var = NULL) {
  ## Split data into a list
  count.list <- list()
  if (!is.null(batch.var)) {
    # TODO: Check if batch variable exists, otherwise break
    
    # Split by batch variable
    if (class(data)[1] == "Seurat") {
      split.object <- SplitObject(data, split.by = batch.var)
      for (ds in 1:length(split.object)) {
        count.list[[ds]] <- Seurat::GetAssayData(split.object[[ds]], slot = "count")
      }
      names(count.list) <- names(split.object)
    } else if (class(data)[1] == "SingleCellExperiment") {
      batch.levels <- unique(colData(data)[,batch.var])
      split.object <- lapply(batch.levels, function(x) data[, colData(data)[,batch.var] == x])
      names(split.object) <- batch.levels
      for (ds in 1:length(split.object)) {
        count.list[[ds]] <- SingleCellExperiment::counts(split.object[[ds]])
      }
      names(count.list) <- names(split.object)
    } else {
      ## TODO: Break with message
    }
  } else {
    # TO DO: Check if data is a list, otherwise break
    for (ds in 1:length(data)) {
      if (class(data[[ds]])[1] == "Seurat") {
        count.list[[ds]] <- Seurat::GetAssayData(data[[ds]], slot = "count")
      } else if (class(data[[ds]])[1] == "SingleCellExperiment") {
        count.list[[ds]] <- SingleCellExperiment::counts(data[[ds]])
      } else if (class(data[[ds]])[1] %in% c("dgCMatrix","dgTMatrix")) {
        count.list[[ds]] <- as(data[[ds]], "dgCMatrix")
      } else if (class(data[[ds]])[1] == "matrix") {
        count.list[[ds]] <- as(data[[ds]], "dgCMatrix")
      } else {
        # TO DO: Verbose and skip
      }
      names(count.list)[ds] <- names(data)[ds]
    }
  }

  ## Return
  return(count.list)
}
