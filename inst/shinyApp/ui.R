library(shiny)
library(shinythemes)

selectDataInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  #UI declaration
  column(
    3,
    h3("Select data source:"),
    #Select data source
    wellPanel(fluidRow(
      column(
        12,
        #Select matrix
        radioButtons("matrixSourceType", "Select matrix:",
                     c("From file" = "fileMatrix",
                       "From workspace" = "varMatrix"))
        ,
        #If the selected data source is a file
        conditionalPanel(
          condition = "input.matrixSourceType == 'fileMatrix'",
          fileInput("matrixFile", "Choose matrix file:",
                    accept = c(
                      "text/csv",
                      "text/comma-separated-values,text/plain",
                      ".csv",".ods",".xls",".xlt")
          )
        ),
        #If the selected data source is a workspace object
        conditionalPanel(
          condition = "input.matrixSourceType == 'varMatrix'",
          selectInput("matrixVar", "Choose matrix var:",
                      ls(envir=.GlobalEnv))
        ),
        fluidRow(column(12,
                        HTML("<br>"))),
        #Select geneset
        radioButtons("genesetSourceType", "Select GeneSet:",
                     c("From file" = "fileGeneset",
                       "From workspace" = "varGeneset"))
        ,
        #If the selected data source is a file
        conditionalPanel(
          condition = "input.genesetSourceType == 'fileGeneset'",
          fileInput("genesetFile", "Choose GeneSet file:",
                    accept = ".gmt")
        ),
        #If the selected data source is a workspace object
        conditionalPanel(
          condition = "input.genesetSourceType == 'varGeneset'",
          selectInput("genesetVar", "Choose GeneSet var:",
                      ls(envir=.GlobalEnv))
        ),
        HTML("<br>"),
        radioButtons("arg", "Control arguments:",
                     c("No" = "no",
                       "Yes" = "yes"))
      )
    ),
    actionButton("button", "Run"))
  )
}

mainDataInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  #UI Definition
  mainPanel(width = 6,
            h2("Generated GSVA data:"),
            textOutput("information"),
            plotOutput("plot"),
            tableOutput("result"),
            uiOutput("download"))
}

argumentsDataInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  #UI Definition
  column(
    3,
    conditionalPanel(
      condition = "input.arg == 'yes'",
      h3("Select arguments:"),
      wellPanel(fluidRow(
        column(
          12,
          selectInput("method", "Choose method:",
                      c("gsva","ssgsea","zscore","plage")),
          selectInput("kcdf", "Choose kcdf:",
                      c("Gaussian","Poisson","none")),
          radioButtons("absRanking", "abs.ranking:",
                       c("False" = FALSE,
                         "True" = TRUE)),
          numericInput("minSz","min.sz:",value = 1),
          numericInput("maxSz","max.sz (Write 0 for infinite):",value = 0),
          numericInput("parallelSz","parallel.sz:",value = 0),
          selectInput("parallelType", "parallel.type:",
                      c("SOCK","MPI","NWS")),
          radioButtons("mxDiff", "mx.diff:",
                       c("True" = TRUE,
                         "False" = FALSE)),
          conditionalPanel(
            condition = "input.method == 'gsva'",
            numericInput("tau1","tau:",value = 1)
          ),
          conditionalPanel(
            condition = "input.method == 'ssgsea'",
            numericInput("tau2","tau:",value = 0.25)
          ),
          conditionalPanel(
            condition = "input.method == 'zscore' || input.method == 'plage'"
          ),
          radioButtons("ssgseaNorm", "ssgsea.norm:",
                       c("True" = TRUE,
                         "False" = FALSE)),
          radioButtons("verbose", "verbose:",
                       c("True" = TRUE,
                         "False" = FALSE))
        )))
    )
  )
}

fluidPage(theme = shinytheme("simplex"),	
	fluidRow(
	  selectDataInput("dataInput"),
	  mainDataInput("mainInput")
	  ,
	  fluidRow(
	    argumentsDataInput("argumentsInput")
	  )
	)
)