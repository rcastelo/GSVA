download_handler <- function(input, output, session) {
  
  output$download <- renderUI({
    tagList(
      downloadButton('downloadData', 'Download'),
      actionButton('closeSave','Save & Close')
    )
  })
  
  #Controls the Download button
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("gsva_es-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      if("matrix" %in% class(generated_gsva)) # if the whole object is a matrix
      {
        dataFrameObs <- as.data.frame(generated_gsva)
        write.csv(dataFrameObs, file)
      }
      else
      {
        if("ExpressionSet" %in% class(generated_gsva)) #If the generated gsva result object is an ExpressionSet
        {
          expressionSetObs <- exprs(generated_gsva)
          dataFrameObs <- as.data.frame(expressionSetObs)
          write.csv(dataFrameObs, file)
        }
        else
        {
          dataFrameObs <- as.data.frame(generated_gsva)
          write.csv(dataFrameObs, file)
        } 
      }
    }
  )
}