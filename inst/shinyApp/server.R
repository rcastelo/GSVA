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
    gsva_information(input,output,session, newY, genes)
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

gsva_information <- function(input, output, session, newY, genes) {
  gsva_es <- NA
  if("matrix" %in% class(generated_gsva))
    gsva_es <- as.data.frame(generated_gsva)
  else if ("ExpressionSet" %in% class(generated_gsva))
    gsva_es <- as.data.frame(exprs(generated_gsva))
  else if ("SummarizedExperiment" %in% class(generated_gsva))
    gsva_es <- as.data.frame(assays(generated_gsva)[[1]])
  else
    stop("Unknown output generated by the call to the 'gsva()' function.")

  #Rendering text1
  output$text1 <- renderUI({
    HTML(paste("<br/>", "\t To see the Empirical Cumulative Distribution Function 
    of a Sample, click on its line in this plot and go
      to the 'Gene.Set' Panel", "<br/>", sep="<br/>"))
  })
  
  # Rendering graph1
  dat.t <- melt(as.data.table(generated_gsva, keep.rownames = "gene.sets"), 
                variable.name = "Sample", id.vars="gene.sets")
  n <- length(levels(dat.t$Sample))
  dd.col <- hcl(h = seq(15, 375, length=n), l = 65, c = 100)[1:n]
  names(dd.col)  <- levels(dat.t$Sample)
  
  
  output$plot <- renderPlotly({
    p <- ggplot(data = dat.t, aes(x=value, color=Sample)) +
      stat_density(geom="line", position = "identity") +
      theme(legend.position = "none") + labs(x="GSVA Scores", y="Density") +
      scale_color_manual("Legend", values = dd.col)
    ggplotly(p, tooltip = "Sample", source = "click1")
  })
  
  # Rendering table
  resultInformation <- matrix(data = c(nrow(generated_gsva),
                                       ncol(generated_gsva)),
                              nrow = 1, ncol = 2)
  colnames(resultInformation) <- c("Nr. of gene sets", "Nr. of samples")
  output$result <- renderTable(resultInformation)
  
  #Rendering text2
  output$text2 <- renderUI({
    title1 <- sample.c()
    h2(tags$b(title1), align ="center")
  })
  
  #Rendering text3
  output$text3 <- renderUI({
    HTML(paste("<br/>", "\t To see the Kernel Density Estimation of genes of 
    any given Gene Set in this Sample,  click on any point in this plot and a
    second plot will appear bellow it", "<br/>", sep="<br/>"))
  })
  
  #Rendering graph2
  eventData1 <- reactive({
    event_data("plotly_click", source = "click1")
  })
  
  sample.c <- reactive({
    req(eventData1())
    ind <- eventData1()$curveNumber+1
    colnames(generated_gsva)[ind]
  })
  
  plot2 <- reactive({
    req(sample.c())
    data <- dat.t[Sample==sample.c()]
    p <- ggplot(data = data, aes(x=value, color=Sample)) +
      stat_ecdf(geom="point") + theme(legend.position = "none") + 
      labs(x="GSVA Scores in selected sample", y="Empirical Cumulative Density") +
      scale_color_manual("Legend", values = dd.col)
    p <- ggplotly(p, source="click2") %>% style(text=data$gene.sets)
  })
  
  output$plot2 <- renderPlotly({
    req(plot2())
    plot2()
  })
  
  # Rendering graph 3
  eventData2 <- reactive({
    event_data("plotly_click", source = "click2")
  })
  
  gene.set <- reactive({
    plot2()$x$data[[1]]$text[eventData2()$pointNumber+1]
  })
  
  output$plot3 <- renderPlotly({
    req(eventData2())
    genes.toplot <- geneIds(genes)[[gene.set()]]
    mt <- match(genes.toplot, rownames(newY))
    x <-  newY[na.omit(mt), sample.c()]
    df <- as.data.frame(x)
    df$x <- as.numeric(df$x)
    df$Gene <- rownames(df)
    df$Sample <- sample.c()
    p1 <- ggplot(data = df, aes(x=x, color = Sample, label = Gene)) +
      stat_density(geom="line", position = "identity") +
      geom_rug() + theme(legend.position = "none") + 
      labs(x="Gene Expressions in selected sample", y="Density") +
      xlim(as.numeric(range(newY))) +
      scale_color_manual("legend", values= dd.col)
    ggplotly(p1, tooltip = c("Gene", "x")) %>% style(hoverinfo="none", traces = 1) %>%
      layout(title = list(text = paste0('<br><sup><i>', gene.set(), '</i></sup>'),
                          font = list(size=15)))
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
