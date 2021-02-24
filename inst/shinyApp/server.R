function(input, output, session) {
  
  errors.gsva <- reactiveVal()
  
  output$errorsGsva <- renderText({
    req(errors.gsva())
    errors.gsva()
  })
  
  observeEvent(input$button, {
    show_modal_spinner(spin="fingerprint", # show the modal window
                       text = "Calculating the GSVA function, this may take a while. 
                       If you launch this session through a console, you can 
                       check your progress there.") 
    
    tryCatch({
      gsva_validation(input,output,session)
      errors.gsva(NULL)
    }, error=function(e){
      errors.gsva(paste0("ERROR=", e))
      return(NULL)
    })
   
    remove_modal_spinner() # remove it when done
    
  })
  

  #Observe the Save & Close button
  observeEvent(input$closeSave, {
    stopApp(generated_gsva) #Stops the app and returns the generated_gsva object
  })
}
