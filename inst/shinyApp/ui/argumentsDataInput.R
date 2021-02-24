argumentsDataInput <- function(id) {
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
          selectInput("method", "Choose method:",
                      c("gsva","ssgsea","zscore","plage")),
          selectInput("kcdf", "Choose kcdf:",
                      c("Gaussian","Poisson","none")),
          radioButtons("absRanking", "abs.ranking:",
                       c("False" = FALSE,
                         "True" = TRUE)),
          numericInput("minSz","min.sz:",value = 1),
          numericInput("maxSz","max.sz (Write 0 for infinite):",value = 0),
          radioButtons("mxDiff", "mx.diff:",
                       c("True" = TRUE,
                         "False" = FALSE)),
          conditionalPanel(
            condition = "input.method == 'gsva'",
            numericInput("tau1","tau:",value = 1)
          ),
          conditionalPanel(
            condition = "input.method == 'ssgsea'",
            numericInput("tau2","tau:",value = 0.25),
            radioButtons("ssgseaNorm", "ssgsea.norm:",
                         c("True" = TRUE,
                           "False" = FALSE)),
          ),
          radioButtons("verbose", "verbose:",
                       c("True" = TRUE,
                         "False" = FALSE))
        )))
    )
  )
}