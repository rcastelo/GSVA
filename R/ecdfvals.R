.ecdfvals_sparse_to_sparse <- function(X, verbose) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  stopifnot(is.logical(verbose)) ## QC
  .Call("ecdfvals_sparse_to_sparse_R", X, Xrsp, verbose)
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

.kcdfvals_sparse_to_sparse <- function(X, Gaussk) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  .Call("kcdfvals_sparse_to_sparse_R", X, Xrsp, Gaussk)
}

.kcdfvals_sparse_to_dense <- function(X, Gaussk) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  Xrsp <- as(X, "RsparseMatrix")
  .Call("kcdfvals_sparse_to_dense_R", X, Xrsp, Gaussk)
}

.order_rankstat <- function(x) {
  stopifnot(is.numeric(x)) ## QC
  .Call("order_rankstat_R", x)
}

.gsva_rnd_walk <- function(gsetIdx, geneRanking, rankStat) {
  stopifnot(is.integer(gsetIdx)) ## QC
  stopifnot(is.integer(geneRanking)) ## QC
  stopifnot(is.integer(rankStat)) ## QC
  .Call("gsva_rnd_walk_R", gsetIdx, geneRanking, rankStat)
}

.gsva_score_genesets <- function(geneSetsRankIdx, geneRanking, rankStat,
                                 maxDiff, absRnk, tau) {
  stopifnot(is.list(geneSetsRankIdx)) ## QC
  stopifnot(length(geneSetsRankIdx) > 0) ## QC
  stopifnot(is.integer(geneSetsRankIdx[[1]])) ## QC
  stopifnot(is.integer(geneRanking)) ## QC
  stopifnot(is.integer(rankStat)) ## QC
  stopifnot(is.logical(maxDiff)) ## QC
  stopifnot(is.logical(absRnk)) ## QC
  stopifnot(is.numeric(tau)) ## QC
  .Call("gsva_score_genesets_R", geneSetsRankIdx, geneRanking, rankStat,
        maxDiff, absRnk, tau)
}

.order_rankstat_sparse_to_dense <- function(X, j) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  .Call("order_rankstat_sparse_to_dense_R", X, j)
}

.order_rankstat_sparse_to_sparse <- function(X, j) {
  stopifnot(is(X, "CsparseMatrix")) ## QC
  .Call("order_rankstat_sparse_to_sparse_R", X, j)
}
