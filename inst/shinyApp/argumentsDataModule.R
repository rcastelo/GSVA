argumentsDataUI <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  #UI Definition
  column(
    width=3,
    conditionalPanel(
      condition = "input.arg == 'yes'",
      h3("Parameters"),
      wellPanel(fluidRow(
        column(
          12,
          selectInput(ns("method"), "Choose method:",
                      c("gsva","ssgsea","zscore","plage")),
          selectInput(ns("kcdf"), "Choose kcdf:",
                      c("Gaussian","Poisson","none")),
          radioButtons(ns("absRanking"), "abs.ranking:",
                       c("False" = FALSE,
                         "True" = TRUE)),
          numericInput(ns("minSz"),"min.sz:",value = 1),
          numericInput(ns("maxSz"),"max.sz (Write 0 for infinite):",value = 0),
          radioButtons(ns("mxDiff"), "mx.diff:",
                       c("True" = TRUE,
                         "False" = FALSE)),
          conditionalPanel(
            condition = "input.method == 'gsva'", ns = ns, 
            numericInput(ns("tau1"),"tau:",value = 1)
          ),
          conditionalPanel(
            condition = "input.method == 'ssgsea'", ns = ns, 
            numericInput(ns("tau2"),"tau:",value = 0.25),
            radioButtons(ns("ssgseaNorm"), "ssgsea.norm:",
                         c("True" = TRUE,
                           "False" = FALSE)))
        )))
    )
  )
}

argumentsDataServer <- function(id){
  moduleServer(id, function(input, output, session){
    varMinsz <-  reactive({
      validate(need(!is.na(input$minSz), "Value 'min.sz' cannot be empty and must be an integer value"))
      input$minSz })
    varMaxsz <- reactive({
      validate(need(!is.na(input$maxSz), "Value 'max.sz' cannot be empty and must be an integer value"))
      ifelse(input$maxSz==0, Inf, input$maxSz) })
    selectedTau <-  reactive({
      if(input$method == "gsva"){
        validate(need(!is.na(input$tau1), "Value 'tau' cannot be empty and must be an integer value"))
        input$tau1
      } else {
        if(input$method == "ssgsea"){
          validate(need(!is.na(input$tau2), "Value 'tau' cannot be empty and must be an integer value"))
          input$tau2
        } else {
          NULL
        }
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