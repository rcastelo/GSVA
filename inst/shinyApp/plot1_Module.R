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
        rv$p <- ggplot(data = rv$dat.t, aes(x=value, color=Sample)) +
          stat_density(geom="line", position = "identity") +
          theme(legend.position = "none") + labs(x="GSVA Scores", y="Density") +
          scale_color_manual("Legend", values = rv$dd.col)
        ggplotly(rv$p, tooltip = "Sample", source = "click1")
      })
    }
  )
}