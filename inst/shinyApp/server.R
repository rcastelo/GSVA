function(input, output, session) {
  
  # CREATE REACTIVE FOR CONSOLE TEXT PROGRESS BAR
  rout <- tempfile("consoleText", fileext = ".txt")
  file.create(rout)
  console.text <- reactiveFileReader(200, session, rout, readLines, warn=F)
  
  # ERRORS MESSAGES
  output$errorsGsva <- renderText({
    req(argInp$varMinsz(), argInp$varMaxsz(), argInp$selectedTau())
    rv$errors.gsva
  })

  # ENABLING 'RUN' BTN
  observe({
    if(!is.null(matrix()) && !is.null(genesets())){
      enable("button")
    } else {
      disable("button")
    }
  })
  
  ### INPUTS ###
  
  # DATA MATRIX
  matrix <- matrixServer("matrix1")
  
  # GENES
  genesets <- geneSetsServer("genes1")
  
  # ARGUMENTS
  argInp <- argumentsDataServer("argumentsInput")

  #### GSVA RESULTS ####
  
  rv <- reactiveValues(gs=NULL, dat.t=NULL, n=NULL, dd.col=NULL, p=NULL, 
                       errors.gsva = NULL, matrix=NULL, genesets=NULL)
  gsva.cancel <- reactiveVal(FALSE)
  
  observeEvent( input$button, {
    rv$gs <- NULL
    rv$dat.t <- NULL
    rv$p <- NULL
    rv$p2 <- NULL
    rv$p3 <- NULL
    rv$errors.gsva = NULL
    rv$matrix <- isolate(matrix())
    rv$genesets <- isolate(genesets())
    gsva.cancel(FALSE)
    modalGSVAUI("modal.text")
    # future() cannot take reactive values, so we must isolate() them
    future({
      sink(rout)
      result <- gsva(isolate(matrix()),
                     isolate(genesets()), 
                     method=isolate(argInp$method()),
                     kcdf=isolate(argInp$kcdf()),
                     abs.ranking=isolate(argInp$absRanking()),
                     min.sz= isolate(argInp$varMinsz()),
                     max.sz=isolate(argInp$varMaxsz()),
                     parallel.sz=1L, ## by now, disable parallelism
                     mx.diff=isolate(argInp$mxDiff()),
                     tau=isolate(argInp$selectedTau()),
                     ssgsea.norm=isolate(argInp$ssgseaNorm()),
                     verbose=TRUE)
      sink()
      write("", file=rout)
      return(result)
    }, seed = TRUE) %...>%
      (function(result){
        rv$gs <- result
        rv$dat.t <- melt(as.data.table(rv$gs, keep.rownames = "gene.sets"),
                         variable.name = "Sample", id.vars="gene.sets")
        rv$n <- length(levels(rv$dat.t$Sample))
        rv$dd.col <- hcl(h = seq(15, 375, length=rv$n), l = 65, c = 100)[1:rv$n]
        names(rv$dd.col)  <- levels(rv$dat.t$Sample)
        write("", file=rout)
        removeModal()
      }) %...!%
      (function(error){
        removeModal()
        write("", file=rout)
        if(gsva.cancel()){
          rv$errors.gsva <- NULL
        } else {
          rv$errors.gsva <- as.character(error)
        }
        
      })
    NULL
  })
  
  # PRINTING CONSOLE.TEXT
  modalGSVAServer("modal.text", console.text, gsva.cancel, rout)
  
  # PLOT1 RENDER
  plot1_Server("plot1", rv)


  # PLOT2 RENDER
  eventData1 <- reactive({
    req(rv$dat.t)
    ind <- event_data("plotly_click", source = "click1")
    ind <- ind$curveNumber+1
  })
  plot2_Server("plot2", eventData1, rv)
  
  
  # PLOT3 RENDER
  eventData2 <- reactive({
    req(rv$p2)
    ind <- event_data("plotly_click", source = "click2")
    ind <- ind$pointNumber+1
  })
  plot3_Server("plot3", eventData2, rv, rv$matrix, rv$genesets)
  
  # DWN BTN
  downloadServer("download", reactive(rv$gs))
  
  # CLOSE BTN
  closeBtnServer("close", reactive(rv$gs))

  
  # TEXT1
  output$text1 <- renderUI({
    req(rv$gs)
    HTML(paste("<br/>", "\t To see the Empirical Cumulative Distribution Function 
    of a Sample, click on its line in this plot and go
      to the 'Gene.Set' Panel", "<br/>", sep="<br/>"))
  })
  
  # TABLE
  output$result <- renderTable({
    req(rv$gs)
    resultInformation <- data.frame("Nr of gene sets" = nrow(rv$gs),
                                    "Nr of samples" = ncol(rv$gs))
    resultInformation
  })
  
  # TEXT2
  output$text2 <- renderUI({
    title1 <- rv$sample.c
    h2(tags$b(title1), align ="center")
  })
  
  # TEXT3
  output$text3 <- renderUI({
    HTML(paste("<br/>", "\t To see the Kernel Density Estimation of genes of 
    any given Gene Set in this Sample,  click on any point in this plot and a
    second plot will appear bellow it", "<br/>", sep="<br/>"))
  })
  
  # SESSION INFO
  output$sessionInfo <- renderPrint({
    sessionInfo()
  })

}
