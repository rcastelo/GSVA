igsva <- function() {

  shinydeps <- c("shiny", "shinydashboard", "ggplot2",
                 "data.table", "plotly", "future",
                 "shinyjs", "shinybusy", "limma", "magrittr")
  maskshinydeps <- shinydeps %in% installed.packages()
  if (any(!maskshinydeps))
    stop(sprintf("Please install the following packages to use the GSVA WebApp:\n\n  %s\n",
                 paste(shinydeps[!maskshinydeps], collapse=", ")))

  appDir <- system.file("shinyApp", package="GSVA")
  if (appDir == "")
    stop("The GSVA Shiny app cannot be found within the package.")


  shiny::runApp(appDir)
}

## legacy code

# runWebApp <- get("runApp", mode="function")
# runWebApp(appDir)


# namespaceImportFrom(self=getNamespace("base"),
#                     ns=getNamespace("shiny"))
# 
# namespaceImportFrom(self=getNamespace("base"),
#                     ns=getNamespace("shinythemes"),
#                     vars="shinytheme")
# 
# namespaceImportFrom(self=getNamespace("base"),
#                     ns=getNamespace("ggplot2"),
#                     vars="ggplot")
# 
# namespaceImportFrom(self=getNamespace("base"),
#                     ns=getNamespace("data.table"))
# 
# namespaceImportFrom(self=getNamespace("base"),
#                     ns=getNamespace("plotly"),
#                     vars=c("ggplotly", "event_data", "style"))
# 
# namespaceImportFrom(self=getNamespace("base"),
#                     ns=getNamespace("future"))
# 
# namespaceImportFrom(self=getNamespace("base"),
#                     ns=getNamespace("shinyjs"))
# 
# namespaceImportFrom(self=getNamespace("base"),
#                     ns=getNamespace("shinybusy"))
# 
# namespaceImportFrom(self=getNamespace("base"),
#                     ns=getNamespace("limma"))