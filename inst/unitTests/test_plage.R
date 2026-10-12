###
### test functions for internal calculations of the PLAGE method
###

test_plageLeadingEigen <- function() {
    message("Running unit tests for the leading eigenpair in PLAGE.")

    ## power iteration is used above 'maxeigen' genes, here 100 instead of
    ## 1000, to run fast; starting from the sum of the columns, which lies
    ## along the vector of ones when the leading eigenvector is a contrast
    ## whose values sum to zero, it converges to the smaller eigenvalue of
    ## the vector of ones, which the Cholesky decomposition rejects, and
    ## eigen() gives the leading one
    k <- 200
    contrast <- rep(c(1, -1), each=k / 2) / sqrt(k)
    ones <- rep(1, k) / sqrt(k)
    gram <- 10 * tcrossprod(contrast) + 3 * tcrossprod(ones) + 0.01 * diag(k)
    e <- GSVA:::.plage_leading_eigen(gram, maxeigen=100)
    checkEqualsNumeric(e$value, 10.01)
    checkEqualsNumeric(abs(sum(e$vector * contrast)), 1)

    ## otherwise, the eigenpair given by power iteration is accepted, and
    ## equals the one of eigen(), up to the sign of the eigenvector
    set.seed(123)
    z <- matrix(rnorm(k * 50), k, 50) + outer(rnorm(k), rnorm(50))
    gram <- tcrossprod(z)
    e <- GSVA:::.plage_leading_eigen(gram, maxeigen=100)
    ref <- eigen(gram, symmetric=TRUE)
    checkEqualsNumeric(e$value, ref$values[1])
    checkEqualsNumeric(abs(sum(e$vector * ref$vectors[, 1])), 1)
}
