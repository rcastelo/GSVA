mainDataInput <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  #UI Definition
  mainPanel(width=6,
            tabsetPanel(type="tabs",
                        tabPanel("Samples",
                                 textOutput("errorsGsva"),
                                 htmlOutput("text1"),
                                 plotlyOutput("plot"),
                                 tableOutput("result"),
                                 uiOutput("download")),
                        tabPanel("Gene Sets",
                                 uiOutput("text2"),
                                 htmlOutput("text3"),
                                 plotlyOutput("plot2"),
                                 plotlyOutput("plot3")),
                        tabPanel("Session Info",
                                 verbatimTextOutput("sessionInfo"))
            )
  )
  
}