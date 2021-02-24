gsva_validation <- function(input, output, session) {

  newY <- if(input$matrixSourceType == "fileMatrix"){
    if(is.null(input$matrixFile))return()
    data.matrix(read.csv(file=input$matrixFile$datapath, row.names = 1L))
  } else {
    get(input$matrixVar)
  }
  
  genes <- 
    if(input$genesetSourceType == "fileGeneset"){
      if(is.null(input$genesetFile))return()
      getGmt(input$genesetFile$datapath)
    } else {
      get(input$genesetVar)
    }
  
  if(input$maxSz == 0) {
    varMaxsz <- Inf
  }else {
    varMaxsz <- input$maxSz
  }
  gsva_generation(input, output, session, newY, genes,varMaxsz)
  gsva_information(input,output,session, newY, genes)
}