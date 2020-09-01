igsva <- function() {

  shinydeps <- c("shiny", "shinythemes", "ggplot2",
                 "data.table", "plotly")
  maskshinydeps <- shinydeps %in% installed.packages()
  if (any(!maskshinydeps))
    stop(sprintf("Please install the following packages to use the GenomiScores WebApp:\n\n  %s\n",
                 paste(shinydeps[!maskshinydeps], collapse=", ")))
  
  namespaceImportFrom(self=getNamespace("base"),
                      ns=getNamespace("shiny"))

  namespaceImportFrom(self=getNamespace("base"),
                      ns=getNamespace("shinythemes"),
                      vars="shinytheme")

  namespaceImportFrom(self=getNamespace("base"),
                      ns=getNamespace("ggplot2"),
                      vars="ggplot")

  namespaceImportFrom(self=getNamespace("base"),
                      ns=getNamespace("data.table"),
                      vars="as.data.table")

  namespaceImportFrom(self=getNamespace("base"),
                      ns=getNamespace("plotly"),
                      vars=c("ggplotly", "event_data", "style"))

  appDir <- system.file("shinyApp", package="GSVA")
  if (appDir == "")
    stop("The GSVA Shiny app cannot be found within the package.")

  runWebApp <- get("runApp", mode="function")
  runWebApp(appDir)
}
