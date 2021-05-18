geneSetsUI <- function(id){
  ns <- NS(id)
  div( id = ns("genesets-input"),
    radioButtons(ns("genesetSourceType"), 
                 label = h5("GENE SETS", style="font-weight: bold"),
                 choices = c("From file" = "fileGeneset",
                   "From workspace" = "varGeneset")),
    conditionalPanel(
      condition = "input.genesetSourceType == 'fileGeneset'", ns = ns,
      fileInput(ns("genesetFile"), "Choose gene sets file:",
                accept = c(".gmt", "text/csv", ".csv"))
    ),
    conditionalPanel(
      condition = "input.genesetSourceType == 'varGeneset'", ns = ns, 
      selectInput(ns("genesetVar"), "Choose gene sets object:",
                  ls(envir=.GlobalEnv))
    )
  )
}

geneSetsServer <- function(id){
  moduleServer( id, function(input, output, session){
    geneSets <- reactive({
      if(input$genesetSourceType == "fileGeneset"){
        if(is.null(input$genesetFile)) return(NULL) #this is in order to disable "run" btn
        ext <- tools::file_ext(input$genesetFile$name)
        genesets <- 
          switch(ext,
                 csv = as.list(read.csv(input$genesetFile$datapath)),
                 gmt = getGmt(input$genesetFile$datapath)
        )
      } else {
        if(is.null(input$genesetVar)) return(NULL)
        genesets <- get(input$genesetVar)
      }
      genesets
    }) %>% bindCache({
      if(input$genesetSourceType == "fileGeneset"){
        input$genesetFile$name
      } else {
        input$genesetVar
      }
    })
  })
}
