closeBtnUI <- function(id){
  ns <- NS(id)
  hidden(actionButton(ns("closeSave"), "SAVE & CLOSE", 
                      icon = icon("window-close"),
                      width = "10vw"))
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