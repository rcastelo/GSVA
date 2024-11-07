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
                  radioButtons(ns("absRanking"), "absRanking:",
                               c("False" = FALSE,
                                 "True" = TRUE)),
                  numericInput(ns("minSz"),"minSize", value = 1),
                  numericInput(ns("maxSz"),"maxSize (Write 0 for infinite)", value = 0),
                  radioButtons(ns("mxDiff"), "maxDiff",
                               c("True" = TRUE,
                                 "False" = FALSE)),
                  numericInput(ns("tau"),"tau", value = 1),
                  radioButtons(ns("ssgseaNorm"), "normalize:",
                               c("True" = TRUE,
                                 "False" = FALSE)),
                  selectInput(ns("checkNA"), "checkNA",
                              c("auto", "yes", "no")),
                  selectInput(ns("use"), "use",
                              c("everything", "all.obs", "na.rm"))
                  )
    )
}

argumentsDataServer <- function(id){
    moduleServer(id, function(input, output, session){
        
        observeEvent(input$method, {
            toggleElement("kcdf", condition = input$method %in% c("gsva"))
            toggleElement("absRanking", condition = input$method %in% "gsva")
            toggleElement("ssgseaNorm", condition = input$method %in% "ssgsea")
            toggleElement("mxDiff", condition = input$method %in% "gsva")
            toggleElement("tau", condition = input$method %in% c("gsva", "ssgsea"))
            toggleElement("checkNA", condition = input$method %in% c("gsva", "ssgsea"))
            toggleElement("use", condition = input$method %in% c("gsva", "ssgsea"))
            
            if(input$method == "gsva"){
                updateNumericInput(inputId = "tau", value = 1)
            } else {
                updateNumericInput(inputId = "tau", value = 0.25)
            }
            
            if(input$method %in% c("zscore", "plage")){
                updateSelectInput(inputId = "kcdf", selected = "Gaussian")
            }
            
        })
        
        varMinsz <-  reactive({
            validate(need(!is.na(input$minSz), "Value 'minSize' cannot be empty and must be an integer value"))
            input$minSz })
        varMaxsz <- reactive({
            validate(need(!is.na(input$maxSz), "Value 'maxSize' cannot be empty and must be an integer value"))
            ifelse(input$maxSz==0, Inf, input$maxSz) })
        selectedTau <-  reactive({
            if(input$method %in% c("gsva", "ssgsea")){
                validate(need(!is.na(input$tau), "Value 'tau' cannot be empty and must be an numeric value"))
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
        checkNA <- reactive({ as.character(input$checkNA) })
        use <- reactive({ as.character(input$use) })
        
        return(list(
            varMinsz = varMinsz,
            varMaxsz = varMaxsz,
            selectedTau = selectedTau,
            method = method,
            kcdf = kcdf,
            absRanking = absRanking,
            mxDiff = mxDiff,
            ssgseaNorm = ssgseaNorm,
            checkNA = checkNA,
            use = use
        ))
    })
}
