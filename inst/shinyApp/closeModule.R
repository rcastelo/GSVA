closeBtnUI <- function(id){
  ns <- NS(id)
  hidden(actionButton(ns("closeSave"), "Save & Close", 
                      icon = icon("window-close"),
                      style = "color: #fff;
                      font-weight: bold;
                      background-color: red;
                      border-color: #fff;"))
}

closeBtnServer <- function(id, gs){
  moduleServer(
    id, 
    function(input, output, session){
      # SAVE & CLOSE BTN
      observe({
        if(is.null(gs())){
          hide("closeSave")
        } else {
          show("closeSave")
        }
      })
      
      observeEvent(input$closeSave, {
        stopApp(gs()) #Stops the app and returns the rv$gs object to the R session
      })
    }
  )
}