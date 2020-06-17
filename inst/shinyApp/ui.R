library(shiny)
library(shinythemes)
library(plotly)

selectDataInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  #UI declaration
  column(
    3,
    h3("Data input:"),
    #Select data source
    wellPanel(fluidRow(
      column(
        12,
        #Select expression data matrix
        radioButtons("matrixSourceType", "Select expression data matrix:",
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
        radioButtons("genesetSourceType", "Select gene sets:",
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
        radioButtons("arg", "Change default settings:",
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
  mainPanel( width = 6,
            tabsetPanel(type="tabs",
                        tabPanel("Samples",
                                 htmlOutput("text1"),
                                 plotlyOutput("plot"),
                                 tableOutput("result"),
                                 uiOutput("download")),
                        tabPanel("Gene Sets",
                                 uiOutput("text2"),
                                 htmlOutput("text3"),
                                 plotlyOutput("plot2"),
                                 plotlyOutput("plot3")),
                        tabPanel("Session Info",
                                 verbatimTextOutput("sessionInfo"))
                                 )
            )
}

argumentsDataInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  #UI Definition
  column(
    3,
    conditionalPanel(
      condition = "input.arg == 'yes'",
      h3("Parameters:"),
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
          radioButtons("mxDiff", "mx.diff:",
                       c("True" = TRUE,
                         "False" = FALSE)),
          conditionalPanel(
            condition = "input.method == 'gsva'",
            numericInput("tau1","tau:",value = 1)
          ),
          conditionalPanel(
            condition = "input.method == 'ssgsea'",
            numericInput("tau2","tau:",value = 0.25),
            radioButtons("ssgseaNorm", "ssgsea.norm:",
                         c("True" = TRUE,
                           "False" = FALSE)),
          ),
          radioButtons("verbose", "verbose:",
                       c("True" = TRUE,
                         "False" = FALSE))
        )))
    )
  )
}

fluidPage(
  theme = shinytheme("spacelab"),	
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  titlePanel(div(h2("GSVA WebApp", align="left"),
             tags$img(src="GSVA.png", align="center", height=75, width=75)),
             windowTitle="GSVA"),
	fluidRow(
	  selectDataInput("dataInput"),
	  mainDataInput("mainInput"),
	  argumentsDataInput("argumentsInput")
	)
)
