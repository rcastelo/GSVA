### LOADING LIBRARIES
library(shiny)
library(plotly)
library(GSVA)
library(GSEABase)
library(limma)
library(ggplot2)
library(data.table)
library(future)
library(promises)
library(shinyjs)
library(shinybusy)
library(shinydashboard)

### there seems to be a problem with the DT package colliding with plotly
### so first check if the package is loaded and if it is, unload it
if("DT" %in% (.packages())) {
  detach("package:DT", unload=TRUE)
}

### setting plan for futures
plan(multisession)

### sourcing all modules
source("argumentsDataModule.R")
source("modalGSVAModule.R")
source("downloadModule.R")
source("plot1_Module.R")
source("plot2_Module.R")
source("plot3_Module.R")
source("matrixModule.R")
source("geneSetsModule.R")
source("argumentsDataModule.R")
source("closeModule.R")

### setting a global in "method" choices for gsva()
methodChoices <- c( "GSVA" = "gsva",
                    "ssGSEA" = "ssgsea",
                    "zscore" = "zscore",
                    "PLAGE" = "plage")