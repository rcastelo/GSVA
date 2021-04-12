fluidPage( 
  theme = shinytheme("spacelab"),
  shinyjs::useShinyjs(),
  add_busy_spinner(spin = "double-bounce", position = "bottom-right", height = "100px", width = "100px"),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  titlePanel(
    fluidRow(
      column(6,
             h2("GSVA Shiny App", align="left")),
      column(6,
             tags$img(src="GSVA.png", align="right", height=75, width=75))
    ), windowTitle="GSVA"),
  
  fluidRow(
    column(
      width=3,
      h3("Data input"),
      #Select data source
      wellPanel(fluidRow(
        column(
          12,
          matrixUI("matrix1"),
          fluidRow(column(12,
                          HTML("<br>"))),
          geneSetsUI("genes1"),
          HTML("<br>"),
          radioButtons("arg", "Change default settings:",
                       c("No" = "no",
                         "Yes" = "yes")),
          actionButton("button", "Run"))
      ))
    ),
    mainPanel(width=6,
              tabsetPanel(type="tabs",
                          tabPanel("Samples",
                                   textOutput("errorsGsva"),
                                   htmlOutput("text1"),
                                   plot1_UI("plot1"),
                                   tableOutput("result"),
                                   downloadUI("download")),
                          tabPanel("Gene Sets",
                                   uiOutput("text2"),
                                   htmlOutput("text3"),
                                   plot2_UI("plot2"),
                                   plot3_UI("plot3")
                                   ),
                          tabPanel("Session Info",
                                   verbatimTextOutput("sessionInfo"))
              )
    ),
    argumentsDataUI("argumentsInput")
    # column(
    #   width=3,
    #   conditionalPanel(
    #     condition = "input.arg == 'yes'",
    #     h3("Parameters"),
    #     wellPanel(fluidRow(
    #       column(
    #         12,
    #         selectInput("method", "Choose method:",
    #                     c("gsva","ssgsea","zscore","plage")),
    #         selectInput("kcdf", "Choose kcdf:",
    #                     c("Gaussian","Poisson","none")),
    #         radioButtons("absRanking", "abs.ranking:",
    #                      c("False" = FALSE,
    #                        "True" = TRUE)),
    #         numericInput("minSz","min.sz:",value = 1),
    #         numericInput("maxSz","max.sz (Write 0 for infinite):",value = 0),
    #         radioButtons("mxDiff", "mx.diff:",
    #                      c("True" = TRUE,
    #                        "False" = FALSE)),
    #         conditionalPanel(
    #           condition = "input.method == 'gsva'", 
    #           numericInput("tau1","tau:",value = 1)
    #         ),
    #         conditionalPanel(
    #           condition = "input.method == 'ssgsea'",
    #           numericInput("tau2","tau:",value = 0.25),
    #           radioButtons("ssgseaNorm", "ssgsea.norm:",
    #                        c("True" = TRUE,
    #                          "False" = FALSE)))
    #       )))
    #   )
    # )
  )
)
