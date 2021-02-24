gsva_generation <- function(input, output, session, newY, genes, varMaxsz) {
  x <- input$method
  selectedTau <- NULL
  switch (x,
          "gsva" = {
            selectedTau <- input$tau1
          },
          "ssgsea" = {
            selectedTau <- input$tau2
          },
          "zscore" = {
            selectedTau <- NULL
          },
          "plage" = {
            selectedTau <- NULL
          }
  )
  # GSVA Generation

  generated_gsva <<- gsva(newY, genes, method=input$method, kcdf=input$kcdf,
                          abs.ranking=as.logical(input$absRanking),
                          min.sz=input$minSz, max.sz=varMaxsz,
                          parallel.sz=1L, ## by now, disable parallelism
                          mx.diff=as.logical(input$mxDiff),
                          tau=selectedTau,
                          ssgsea.norm=as.logical(input$ssgseaNorm),
                          verbose=as.logical(input$verbose))
}