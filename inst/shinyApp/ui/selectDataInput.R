selectDataInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  #UI declaration
  column(
    width=3,
    h3("Data input"),
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