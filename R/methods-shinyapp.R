#'
#' @importFrom shiny HTML actionButton animationOptions checkboxGroupInput column div downloadHandler downloadLink eventReactive fileInput fluidPage fluidRow h2 h3 h4 headerPanel htmlOutput mainPanel need numericInput NS observe observeEvent p plotOutput reactiveValues renderPlot renderUI selectInput shinyApp sliderInput stopApp tabPanel tabsetPanel textOutput uiOutput updateSelectInput validate wellPanel withProgress conditionalPanel reactive outputOptions tableOutput tags radioButtons downloadButton
#' @importFrom shinythemes shinytheme
#' @importFrom utils head
#' @importFrom geneplotter multidensity
#' @importFrom stats median
#' @importFrom graphics plot
#' @export
#'

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

gsva_validation <- function(input, output, session) {
  success <- FALSE #Variable to control if the GSVA variables are assigned correctly
  if(input$matrixSourceType == "fileMatrix")
  {
    if (is.null(input$matrixFile))
    {
      paste("No matrix file selected!")
      success <- FALSE
    }
    else
    {
      #Matrix file selected
      if(input$genesetSourceType == "fileGeneset")
      {
        if (is.null(input$genesetFile))
        {
          paste("No geneSet file selected!")
          success <- FALSE
        }
        else
        {
          #User selects matrix file and geneSet file
          inFile <- input$matrixFile
          newY <- as.matrix(read.csv(inFile$datapath, header = TRUE, sep = ",")) #Reading file as matrix
          rownames(newY) <- newY[,1] #Taking the first column as rownames
          newY <- newY[,-1] #Deleting the first column
          inGenesetFile <- input$genesetFile
          genes <- getGmt(inGenesetFile$datapath)
          if(input$maxSz == 0) {
            varMaxsz <- Inf
          }else {
            varMaxsz <- input$maxSz
          }
          success <- TRUE
        }
      }
      else
      {
        #User selects matrix file and geneset var
        inFile <- input$matrixFile
        newY <- as.matrix(read.csv(inFile$datapath, header = TRUE, sep = ",")) #Reading file as matrix
        rownames(newY) <- newY[,1] #Taking the first column as rownames
        newY <- newY[,-1] #Deleting the first column
        assign("genes",get(input$genesetVar))
        if(input$maxSz == 0) {
          varMaxsz <- Inf
        }else {
          varMaxsz <- input$maxSz
        }
        success <- TRUE
      }
    }
  }
  else
  {
    #User selects matrix varand geneset file
    if(input$genesetSourceType == "fileGeneset")
    {
      if (is.null(input$genesetFile))
      {
        paste("No geneSet file selected!")
        success <- FALSE
      }
      else
      {
        assign("newY",get(input$matrixVar))
        inGenesetFile <- input$genesetFile
        genes <- getGmt(inGenesetFile$datapath)
        if(input$maxSz == 0) {
          varMaxsz <- Inf
        }else {
          varMaxsz <- input$maxSz
        }
        success <- TRUE
      }
    }
    else
    {
      #User selects matrix var selected and geneset var
      assign("newY",get(input$matrixVar))
      assign("genes",get(input$genesetVar))
      if(input$maxSz == 0) {
        varMaxsz <- Inf
      }else {
        varMaxsz <- input$maxSz
      }
      success <- TRUE
    }
  }
  if(success==TRUE)
  {
    gsva_generation(input, output, session, newY, genes,varMaxsz)
    gsva_information(input,output,session)
  }
}

gsva_generation <- function(input, output, session, newY, genes,varMaxsz) {
  x <- input$method
  selectedTau <- NULL
  switch (x,
          "gsva" = {
            selectedTau <- input$tau1
          },
          "ssgsea" = {
            selectedTau <- input$tau2
          },
          "zscore" = {
            selectedTau <- NULL
          },
          "plage" = {
            selectedTau <- NULL
          }
  )
  #GSVA Generation
  withProgress(message = 'Runing GSVA', value = 0, {
    incProgress(1, detail = "This may take a while...")
    generated_gsva <<- gsva(newY, genes, method=input$method, kcdf=input$kcdf, abs.ranking=as.logical(input$absRanking),
                            min.sz=input$minSz, max.sz=varMaxsz, parallel.sz=input$parallelSz, parallel.type=input$parallelType,
                            mx.diff=as.logical(input$mxDiff), tau=selectedTau, ssgsea.norm=as.logical(input$ssgseaNorm),
                            verbose=as.logical(input$verbose))
  })
}

gsva_information <- function(input, output, session) {
  if(class(generated_gsva) == "matrix")
  {
    resultInformation <- matrix(data = c(input$matrixVar,input$genesetVar,ncol(generated_gsva),nrow(generated_gsva)), nrow = 1, ncol = 4)
    colnames(resultInformation) <- c("Matrix used","GeneSet used", "Col num", "Row num")
    output$result <- renderTable(resultInformation)
    output$plot <- renderPlot(multidensity(as.list(as.data.frame(generated_gsva)), legend=NA, las=1, xlab=sprintf("%s scores", input$method), main="", lwd=2))
    tagList(
      downloadButton('downloadData', 'Download'),
      actionButton('closeSave','Save & Close')
    )
  }
  else
  {

    resultInformation <- matrix(data = c(input$matrixVar,input$genesetVar,ncol(generated_gsva),nrow(generated_gsva)), nrow = 1, ncol = 4)
    colnames(resultInformation) <- c("Matrix used","GeneSet used", "Col num", "Row num")
    output$result <- renderTable(resultInformation)
    if(class(generated_gsva) == "ExpressionSet") #If the generated gsva is an ExpressionSet
    {
      expressionSetObs <- exprs(generated_gsva)
      output$plot <- renderPlot(multidensity(as.list(as.data.frame(expressionSetObs)), legend=NA, las=1, xlab=sprintf("%s scores", input$method), main="", lwd=2)) 
    }
    else
    {
      output$plot <- renderPlot(multidensity(as.list(as.data.frame(generated_gsva)), legend=NA, las=1, xlab=sprintf("%s scores", input$method), main="", lwd=2))
    }
    tagList(
      downloadButton('downloadData', 'Download'),
      actionButton('closeSave','Save & Close')
    )
  }
}

download_handler <- function(input, output, session) {
  #Controls the Download button
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("gsva_es-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      if(class(generated_gsva) == "matrix") #If the whole object is a matrix
      {
        dataFrameObs <- as.data.frame(generated_gsva)
        write.csv(dataFrameObs, file)
      }
      else
      {
        if(class(generated_gsva) == "ExpressionSet") #If the generated gsva result value is an ExpressionSet
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

igsva <- function() {
  app <- list(ui = NULL, server = NULL)
  app$ui <- fluidPage(theme = shinytheme("simplex"),
                      fluidRow(
                        selectDataInput("dataInput"),
                        mainDataInput("mainInput")
                        ,
                        fluidRow(
                          argumentsDataInput("argumentsInput")
                        )
                      )
  )
  
  app$server <- function(input, output, session) {
    v <- reactiveValues(action = FALSE)
    
    observeEvent(input$button, {
      v$action <- input$button
    })
    
    output$download <- renderUI({
      if(v$action)
      {
        #Isolates the Run event, that allows the program to run the generation only if the user clicks the button.
        isolate({
          gsva_validation(input,output,session)
        })
      }
    })
    download_handler(input,output,session)
    
    #Observe the Save & Close button
    observeEvent(input$closeSave, {
      stopApp(generated_gsva) #Stops the app and returns the generated_gsva object
    })
  }
  runApp(app)
}
