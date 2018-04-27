igsva <- function() {
  appDir <- system.file("shinyApp", package="GSVA")
  if (appDir == "")
    stop("The GSVA Shiny app cannot be found within the package.")

  runApp(appDir)
}
