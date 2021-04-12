geneSetsUI <- function(id){
  ns <- NS(id)
  tagList(
    radioButtons(ns("genesetSourceType"), "Select gene sets:",
                 c("From file" = "fileGeneset",
                   "From workspace" = "varGeneset")),
    conditionalPanel(
      condition = "input.genesetSourceType == 'fileGeneset'", ns = ns,
      fileInput(ns("genesetFile"), "Choose GeneSet file:",
                accept = ".gmt")
    ),
    conditionalPanel(
      condition = "input.genesetSourceType == 'varGeneset'", ns = ns, 
      selectInput(ns("genesetVar"), "Choose GeneSet var:",
                  ls(envir=.GlobalEnv))
    )
  )
}

geneSetsServer <- function(id){
  moduleServer( id, function(input, output, session){
    geneSets <- reactive({
      if(input$genesetSourceType == "fileGeneset"){
        if(is.null(input$genesetFile)) return(NULL)
        genesets <- getGmt(input$genesetFile$datapath)
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