library(GSVA)
library(shiny)
library(shinythemes)
library(GSEABase)
library(GSVAdata)
library(limma)
library(ggplot2)
library(data.table)
library(plotly)

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
          ## numericInput("parallelSz","parallel.sz:",value = 0),
          ## selectInput("parallelType", "parallel.type:",
          ##             c("SOCK","MPI","NWS")),
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
    #User selects matrix var and geneset file
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
  # GSVA Generation
  withProgress(message = 'Runing GSVA', value = 0, {
    incProgress(1, detail = "This may take a while...")
    generated_gsva <<- gsva(newY, genes, method=input$method, kcdf=input$kcdf,
                            abs.ranking=as.logical(input$absRanking),
                            min.sz=input$minSz, max.sz=varMaxsz,
                            parallel.sz=1L, ## by now, disable parallelism
                            mx.diff=as.logical(input$mxDiff),
                            tau=selectedTau,
                            ssgsea.norm=as.logical(input$ssgseaNorm),
                            verbose=as.logical(input$verbose))
  })

}

gsva_information <- function(input, output, session) {
  gsva_es <- NA
  if("matrix" %in% class(generated_gsva))
    gsva_es <- as.data.frame(generated_gsva)
  else if ("ExpressionSet" %in% class(generated_gsva))
    gsva_es <- as.data.frame(exprs(generated_gsva))
  else if ("SummarizedExperiment" %in% class(generated_gsva))
    gsva_es <- as.data.frame(assays(generated_gsva)[[1]])
  else
    stop("Unknown output generated by the call to the 'gsva()' function.")
  
  # Rendering table
  resultInformation <- matrix(data = c(nrow(generated_gsva),
                                       ncol(generated_gsva)),
                                       nrow = 1, ncol = 2)
  colnames(resultInformation) <- c("Nr. of gene sets", "Nr. of samples")
  output$result <- renderTable(resultInformation)
  
  # Rendering graph
  dat.t <- melt(as.data.table(generated_gsva), variable.name = "Sample")
  output$plot <- renderPlotly({
    p <- ggplot(data = dat.t, aes(x=value, color=Sample)) +
      stat_density(geom="line", position = "identity") +
      theme(legend.position = "none") + labs(x="GSVA Scores", y="Density")
    ggplotly(p, tooltip = "Sample")
  })
  
  # Rendering Session Info
  output$sessionInfo <- renderPrint({
    sessionInfo()
  })
  
  tagList(
    downloadButton('downloadData', 'Download'),
    actionButton('closeSave','Save & Close')
  )
}

download_handler <- function(input, output, session) {
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

function(input, output, session) {
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
