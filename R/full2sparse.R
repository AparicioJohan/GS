

#' full2sparse
#'
#' @param mat
#'
#' @return
#' @export
#'
#' @examples
"full2sparse" <- function(mat) {
  if (is.null(colnames(mat))) {
    colnames(mat) <- 1:ncol(mat)
  }
  regions <- as.numeric(colnames(mat))
  sparse <- data.table::as.data.table(expand.grid(regions, regions))
  up.triangle <- c(upper.tri(mat, diag = TRUE))
  sparse <- cbind(sparse, c(mat), up.triangle)
  colnames(sparse) <- c("region1", "region2", "IF", "up.triangle")
  sparse <- subset(sparse, up.triangle == TRUE & IF != 0)
  sparse <- sparse[, 1:3, with = FALSE]
  sparse <- sparse[,c(2,1,3)]
  return(sparse)
}
