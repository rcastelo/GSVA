#'
#' @importFrom shiny HTML actionButton animationOptions checkboxGroupInput column div downloadHandler downloadLink eventReactive fileInput fluidPage fluidRow h2 h3 h4 headerPanel htmlOutput mainPanel need numericInput NS observe observeEvent p plotOutput reactiveValues renderPlot renderUI selectInput shinyApp sliderInput stopApp tabPanel tabsetPanel textOutput uiOutput updateSelectInput validate wellPanel withProgress conditionalPanel reactive outputOptions tableOutput tags radioButtons downloadButton
#' @importFrom shinythemes shinytheme
#' @importFrom utils head
#' @importFrom geneplotter multidensity
#' @importFrom stats median
#' @importFrom graphics plot
#' @export
#'

igsva <- function() {
  runApp("GSVA/R/app")
}