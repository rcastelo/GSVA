closeBtnUI <- function(id){
  ns <- NS(id)
  hidden(actionButton(ns("closeSave"), "Save & Close", 
                      icon = icon("window-close"),
                      width = "10vw",
                      style = "color: #fff; 
                      background-color: red;
                      font-weight: bold;
                      border-color: #fff;
                      padding: 5px 5px 5px 5px;
                      margin: 6px 5px 6px 15px;"
                      ))
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