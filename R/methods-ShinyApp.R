
#' @title Gene Set Variation Analysis
#' 
#' @description Starts an interactive GSVA shiny web app.
#' 
#' GSVA assesses the relative enrichment of gene sets across samples using a
#' non-parametric approach.  Conceptually, GSVA transforms a p-gene by n-sample
#' gene expression matrix into a g-geneset by n-sample pathway enrichment
#' matrix. This facilitates many forms of statistical analysis in the 'space'
#' of pathways rather than genes, providing a higher level of interpretability.
#' 
#' The `igsva()` function starts an interactive shiny web app that allows
#' the user to configure the arguments of the [`gsva()`] function and
#' runs it on the computer. Please see the manual page of the
#' [`gsva()`] function for a description of the arguments and their
#' default and alternative values.
#' 
#' The input data may be loaded from the users workspace or by selecting a CSV
#' file for the expression data, and a GMT file for the gene sets data.
#' 
#' @usage igsva()
#' 
#' @return A gene-set by sample matrix of GSVA enrichment scores after pressing
#' the button 'Save & Close'. This result can be also downloaded as a CSV file
#' with the 'Download' button.
#' 
#' @author J. Fernández and R. Castelo
#' 
#' @seealso [`gsva()`]
#'
#' @aliases igsva
#' 
#' @name igsva
#' 
#' @rdname igsva
#' 
#' @references Hänzelmann, S., Castelo, R. and Guinney, J. GSVA: Gene set
#' variation analysis for microarray and RNA-Seq data.
#' *BMC Bioinformatics*, 14:7, 2013.
#' 
#' @keywords GSVA shiny
#' @examples
#' \dontrun{
#' res <- igsva() ## this will open your browser with the GSVA shiny web app
#' }
#'
#' @export
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
