dashboardPage(
  title = "GSVA Shiny Application",
  dashboardHeader(
    tags$li(class = "dropdown",
            tags$style(".main-header {max-height: 75px}"),
            tags$style(".main-header .logo {height: 75px}"),
            tags$div("GSVA Shiny Application", style = "font-size: 30px; 
                     color: white; font-weight: bold;")
    ),
    title = tags$img(src="GSVA.png", height=75, width=75)
  ),
  dashboardSidebar(
    tags$style(".left-side, .main-sidebar {padding-top: 75px}"),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    h3("Data input"),
    #Select data source
    matrixUI("matrix1"),
    br(),
    geneSetsUI("genes1"),
    br(),
    radioButtons("arg", "Change default settings:",
                 c("No" = "no",
                   "Yes" = "yes")),
    actionButton("button", "Run"),
    br(),
    downloadUI("download"),
    closeBtnUI("close")
  ),
  dashboardBody(
    shinyjs::useShinyjs(),
    add_busy_spinner(spin = "double-bounce", position = "bottom-right",
                     height = "100px", width = "100px"),
    box(
      width = 9,
      tabsetPanel(id = "Panels", type="tabs",
                  tabPanel("Samples",
                           textOutput("errorsGsva"),
                           htmlOutput("text1"),
                           plot1_UI("plot1"),
                           tableOutput("result")
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
