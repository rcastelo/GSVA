plot2_UI <- function(id){
  ns <- NS(id)
  plotlyOutput(ns("plot2"))
}

plot2_Server <- function(id, eventData1, rv){
  moduleServer(
    id, 
    function(input, output, session){
      
      output$plot2 <- renderPlotly({
        req(eventData1())
        
        # in order to print the name of the method (and not the 
        # selected value from the method) on the 'x' label, this
        # name is retrieved from the list 'methodChoices' declared
        # in 'global.R
        method <- names(methodChoices)[methodChoices == rv$method]
        
        rv$sample.c <- colnames(rv$gs)[eventData1()]
        data <- rv$dat.t[Sample==rv$sample.c]
        p <- ggplot(data = data, aes(x=value, color=Sample)) +
          stat_ecdf(geom="point") + theme(legend.position = "none") +
          labs(x=paste0(method, " scores in selected sample"), y="Empirical cumulative distribution") +
          scale_color_manual("Legend", values = rv$dd.col)
        rv$p2 <- ggplotly(p, source="click2") %>% style(text=data$gene.sets)
        rv$p2
      })
    }
  )
}
