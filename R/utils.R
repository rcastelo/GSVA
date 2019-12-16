.filterFeatures <- function(expr, method) {

  ## filter out genes with constant expression values
  sdGenes <- apply(expr, 1, sd)
  if (any(sdGenes == 0) || any(is.na(sdGenes))) {
    warning(sum(sdGenes == 0 | is.na(sdGenes)),
            " genes with constant expression values throuhgout the samples.")
    if (method != "ssgsea") {
      warning("Since argument method!=\"ssgsea\", genes with constant expression values are discarded.")
      expr <- expr[sdGenes > 0 & !is.na(sdGenes), ]
    }
  } 

  expr
}
