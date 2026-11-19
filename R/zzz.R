.onLoad <- function(libname, pkgname) {
  options(Matrix.warnDeprecatedCoerce = 2)
  ## the R process that loads GSVA, whose forked workers inherit it, see
  ## .release_block_mem()
  gsva_global$pid <- Sys.getpid()
}
