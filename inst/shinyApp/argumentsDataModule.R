argumentsDataUI <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  #UI Definition
  conditionalPanel(
    condition = "input.arg == 'yes'",
    fluidRow(
      column(
        width = 12,
        align = "center",
        h4("PARAMETERS", style="font-weight: bold")
      )
    ),
    wellPanel( id= "args-well",
      selectInput(ns("method"), "Method",
                  choices = methodChoices),
      selectInput(ns("kcdf"), "kcdf",
                  c("Gaussian","Poisson","none")),
      radioButtons(ns("absRanking"), "abs.ranking:",
                   c("False" = FALSE,
                     "True" = TRUE)),
      numericInput(ns("minSz"),"min.sz", value = 1),
      numericInput(ns("maxSz"),"max.sz (Write 0 for infinite)", value = 0),
      radioButtons(ns("mxDiff"), "mx.diff",
                   c("True" = TRUE,
                     "False" = FALSE)),
      numericInput(ns("tau"),"tau", value = 1),
      radioButtons(ns("ssgseaNorm"), "ssgsea.norm:",
                   c("True" = TRUE,
                     "False" = FALSE))
    )
  )
}

argumentsDataServer <- function(id){
  moduleServer(id, function(input, output, session){
    
    observeEvent(input$method, {
      toggleElement("kcdf", condition = input$method %in% c("gsva", "ssgsea"))
      toggleElement("absRanking", condition = input$method %in% "gsva")
      toggleElement("ssgseaNorm", condition = input$method %in% "ssgsea")
      toggleElement("mxDiff", condition = input$method %in% "gsva")
      toggleElement("tau", condition = input$method %in% c("gsva", "ssgsea"))
      
      if(input$method == "gsva"){
        updateNumericInput(inputId = "tau", value = 1)
      } else {
        updateNumericInput(inputId = "tau", value = 0.25)
      }
      
      if(input$method %in% c("zscore", "plage")){
        updateSelectInput(inputId = "kcdf", selected = "Gaussian")
      }
      
    })
    
    #"absRanking", "ssgseaNorm", "mxDiff", "tau"
    varMinsz <-  reactive({
      validate(need(!is.na(input$minSz), "Value 'min.sz' cannot be empty and must be an integer value"))
      input$minSz })
    varMaxsz <- reactive({
      validate(need(!is.na(input$maxSz), "Value 'max.sz' cannot be empty and must be an integer value"))
      ifelse(input$maxSz==0, Inf, input$maxSz) })
    selectedTau <-  reactive({
      if(input$method %in% c("gsva", "ssgsea")){
        validate(need(!is.na(input$tau), "Value 'tau' cannot be empty and must be an integer value"))
        input$tau
      } else {
        NULL
      }
    })
    method <-  reactive({ input$method })
    kcdf <-  reactive({ input$kcdf })
    absRanking <-   reactive({ as.logical(input$absRanking) })
    mxDiff <-   reactive({ as.logical(input$mxDiff) })
    ssgseaNorm <-  reactive({ as.logical(input$ssgseaNorm) })
    
    return(list(
      varMinsz = varMinsz,
      varMaxsz = varMaxsz,
      selectedTau = selectedTau,
      method = method,
      kcdf = kcdf,
      absRanking = absRanking,
      mxDiff = mxDiff,
      ssgseaNorm = ssgseaNorm
    ))
  })
}