function(input, output, session) {
  
  # CREATE REACTIVE FOR CONSOLE TEXT PROGRESS BAR
  rout <- tempfile("consoleText", fileext = ".txt")
  file.create(rout)
  console.text <- reactiveFileReader(200, session, rout, readLines, warn=F)
  
  
  ##################### INPUTS  ##################### 
  
  # DATA MATRIX
  matrix <- matrixServer("matrix1")
  
  # GENES
  genesets <- geneSetsServer("genes1")
  
  # ARGUMENTS
  argInp <- argumentsDataServer("argumentsInput")
  

  ##################### GSVA RESULTS  ##################### 
  
  ## REACTIVE VALUES
  rv <- reactiveValues(gs=NULL, dat.t=NULL, n=NULL, dd.col=NULL, p=NULL, 
                       p2=NULL, p3=NULL, errors.gsva = NULL, sample.c = NULL,
                       method=NULL)
  gsva.cancel <- reactiveVal(FALSE)
  
  ## GSVA RESULT
  observeEvent( input$button, {
    
    ## This js is in order to reset the event_data from the plotlys,
    ## so every time the .user hits the 'run' button, plotlys get back to null
    runjs("Shiny.setInputValue('plotly_click-click1', null);")
    runjs("Shiny.setInputValue('plotly_click-click2', null);")
    
    ## here we reset all the reactiveValues to NULL
    rv$gs <- NULL
    rv$dat.t <- NULL
    rv$p <- NULL
    rv$p2 <- NULL
    rv$p3 <- NULL
    rv$sample.c <- NULL
    rv$errors.gsva <- NULL
    rv$method <- argInp$method()
    
    ## this is a flag for the future. Futures cannot be canceled or
    ## terminated in a strict way, so when they get interrupted they
    ## throw an error that is not related to gsva(). When future is 
    ## interrupted, the flag goes TRUE in order to make the errors
    ## message print NULL
    gsva.cancel(FALSE)
    
    modalGSVAUI("modal.text")
    
    ## future() cannot take reactive values, so we must isolate() them
    future({
      ## sink() will redirect all console cats and prints to a
      ## text file that the main session will be reading in order
      ## to print the progress bar from bplaply()
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
      ## when gsva() ends, we reset the console text file to empty
      write("", file=rout)
      return(result)
    }, seed = TRUE) %...>%
      (function(result){
        ## the future's result will be the gsva() result, and we save it
        ## and transform it in reactiveValues(). In order to make the future
        ## not block the app at an inner-session level, we save the results in
        ## reactiveValues() and then at the end of the observeEvent() we return NULL
        ## in order to make the plots.
        ## https://github.com/rstudio/promises/issues/23#issuecomment-386687705
        rv$gs <- result
        rv$dat.t <- melt(as.data.table(rv$gs, keep.rownames = "gene.sets"),
                         variable.name = "Sample", id.vars="gene.sets")
        rv$n <- length(levels(rv$dat.t$Sample))
        rv$dd.col <- hcl(h = seq(15, 375, length=rv$n), l = 65, c = 100)[1:rv$n]
        names(rv$dd.col)  <- levels(rv$dat.t$Sample)
        
        ## finally, we leave the console.text file empty again and
        ## remove the modal
        write("", file=rout)
        removeModal()
      }) %...!%
      (function(error){
        ## there can be two ways to get an error here: 
        ## 1. gsva() fails, which is an ok error and should be returnet to user
        ## 2. User interrupts the future, which shouldn't be printed, that's
        ## why I use a flag to identify if error comes from pressing "Cancel" btn
        ## on the modal
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
  
  
  ##################### OUTPUTS ##################
  
  # PLOT1 RENDER
  plot1_Server("plot1", rv)

  # PLOT2 RENDER
  eventData1 <- reactive({
    if(is.null(rv$p))return(NULL)
    ind <- event_data("plotly_click", source = "click1")
    ind <- ind$curveNumber+1
  })
  plot2_Server("plot2", eventData1, rv)

  # PLOT3 RENDER
  
  ## Whenever the user clicks on the first plot, the third one resets
  observeEvent(eventData1(), {
    runjs("Shiny.setInputValue('plotly_click-click2', null);")
  })
  
  eventData2 <- reactive({
    req(rv$p2)
    ind <- event_data("plotly_click", source = "click2")
    ind <- ind$pointNumber+1
  })
  plot3_Server("plot3", eventData2, rv, matrix, genesets)

  # ERRORS MESSAGES
  output$errorsGsva <- renderText({
    req(argInp$varMinsz(), argInp$varMaxsz(), argInp$selectedTau())
    rv$errors.gsva
  })
  
  # SESSION INFO
  output$sessionInfo <- renderPrint({
    sessionInfo()
  })
  
  
  ##################### UI SETUPS #####################
  
  ## ENABLING 'RUN' BTN
  observe({
    if(!is.null(matrix()) && !is.null(genesets())){
      enable("button")
    } else {
      disable("button")
    }
  })
  
  ## HIDE 'GeneSets' PANEL WHILE THERE IS NO CLICK EVENT ON THE FIRST PLOT
  observe({
    if( length(eventData1()) == 0){
      hideTab(inputId = "Panels", target = "GeneSets")
    } else {
      showTab(inputId = "Panels", target = "GeneSets", select = TRUE)
    }
  })
  
  # DNLD BTN
  downloadServer("download", reactive(rv$gs))
  
  # CLOSE BTN
  closeBtnServer("close", reactive(rv$gs))
  
  
  # TEXT1
  output$text1 <- renderUI({
    req(rv$gs)
    tagList(
      br(),
      div("Non-parametric kernel density estimation of sample
          profiles of GSVA enrichment scores. Clicking on the
          line of a sample will display the empirical cumulative
          distribution of GSVA scores for that sample on the
          'GeneSets' tab", style="text-align: center;")
    )
  })
  
  # TABLE
  output$result <- renderTable({
    req(rv$gs)
    resultInformation <- data.frame("Nr. of gene sets" = nrow(rv$gs),
                                    "Nr. of samples" = ncol(rv$gs),
                                    check.names=FALSE)
    resultInformation
  }, bordered = TRUE)
  
  # TEXT2
  output$text2 <- renderUI({
    title1 <- rv$sample.c
    h2(tags$b(title1), align ="center")
  })
  
  # TEXT3
  output$text3 <- renderUI({
    tagList(
      br(),
      div("Empirical cumulative distribution of GSVA scores, where each
          point is a gene set. Clicking on a gene set will display below
          the individual gene expression values of its constituent genes
          and the non-parametric kernel density estimation of their
          distribution", style = "text-align: center;")
    )
  })
  
}
