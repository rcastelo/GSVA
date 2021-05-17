dashboardPage(
  title = "GSVA Shiny Application",
  
  dashboardHeader(
    tags$li(class = "dropdown",
            tags$div(id = "app_title", "GSVA Shiny Application")
    ),
    title = tags$img(src="GSVA.png", height=75, width=75)
  ),
  
  dashboardSidebar(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    div(h3("DATA INPUT", style="font-weight: bold"), align = "center"),
    br(),
    matrixUI("matrix1"),
    br(),
    geneSetsUI("genes1"),
    br(),
    radioButtons(inputId = "arg",
                 label = h5("CHANGE DEFAULT SETTINGS?", style="font-weight: bold"),
                 c("No" = "no",
                   "Yes" = "yes")),
    br(),
    fluidRow(
      column(
        width = 12, align = "left",
        actionButton("button", "RUN", class = "run-btn", icon = icon("play-circle"),
                     width = "10vw"),
        downloadUI("download"),
        closeBtnUI("close")
      )
    )
  ),
  
  dashboardBody(
    shinyjs::useShinyjs(),
    add_busy_spinner(spin = "cube-grid", position = "bottom-right",
                     height = "100px", width = "100px"),
    fluidRow(
      box(
        width = 9,
        tabsetPanel(id = "Panels", type="tabs",
                    tabPanel("Samples",
                             textOutput("errorsGsva"),
                             htmlOutput("text1"),
                             plot1_UI("plot1"),
                             br(),
                             fluidRow(
                               column(
                                 width = 12,
                                 align = "center",
                                 tableOutput("result")
                               )
                             )
                    ),
                    tabPanel("GeneSets",
                             uiOutput("text2"),
                             htmlOutput("text3"),
                             plot2_UI("plot2"),
                             plot3_UI("plot3")
                    ),
                    tabPanel("Session Info",
                             verbatimTextOutput("sessionInfo"))
        )
      ),
      box(
        width = 3,
        argumentsDataUI("argumentsInput")
      )
    )
  )
  
)
