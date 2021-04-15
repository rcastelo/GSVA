modalGSVAUI <- function(id){
  ns <- NS(id)
  showModal(
    modalDialog(
      size="l",
      title = "GSVA calculations",
      div(id="install.text", "Calculating, please wait... "),
      div(p("\n")),
      verbatimTextOutput(ns("text")),
      footer = actionButton(ns("cancel"), "Cancel"))
  )
}

modalGSVAServer <- function(id, console.text, gsva.cancel, rout){
  moduleServer(
    id,
    function(input, output, session){

      output$text <- renderText({
        req(console.text())
        max <- length(console.text())
        if(max>1){
          paste(console.text()[1], console.text()[max], sep= "\n")
        } else {
          console.text()
        }
      })
      
      observeEvent(input$cancel, {
        removeModal()
        gsva.cancel(TRUE)
        write("", file=rout)
        # changing plan() is the only and recommended way to safely 'interrupt' a future process:
        # https://github.com/HenrikBengtsson/future/issues/93
        plan(sequential)
        plan(multisession)
      })
    }
  )
}