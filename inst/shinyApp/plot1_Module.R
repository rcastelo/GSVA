plot1_UI <- function(id){
  ns <- NS(id)
  plotlyOutput(ns("plot"))
}

plot1_Server <- function(id, rv){
  moduleServer(
    id,
    function(input, output, session){
      
      output$plot <- renderPlotly({
        req(rv$dat.t)
        # in order to print the name of the method (and not the 
        # selected value from the method) on the 'x' label, this
        # name is retrieved from the list 'methodChoices' declared
        # in 'global.R
        method <- names(methodChoices)[methodChoices == rv$method]
        rv$p <- ggplot(data = rv$dat.t, aes(x=value, color=Sample)) +
          stat_density(geom="line", position = "identity") +
          theme(legend.position = "none") + labs(x=paste0(method, " Scores"), y="Density") +
          scale_color_manual("Legend", values = rv$dd.col)
        ggplotly(rv$p, tooltip = "Sample", source = "click1")
      })
    }
  )
}