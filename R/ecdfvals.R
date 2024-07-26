.ecdfvals_sparse_to_sparse <- function(X) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  .Call("ecdfvals_sparse_to_sparse_R", X, Xrsp)
}

.ecdfvals_sparse_to_dense <- function(X) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  .Call("ecdfvals_sparse_to_dense_R", X, Xrsp)
}

.ecdfvals_dense_to_dense <- function(X) {
  stopifnot(is.matrix(X)) ## QC
  .Call("ecdfvals_dense_to_dense_R", X)
}
