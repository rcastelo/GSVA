plot3_UI <- function(id){
  ns <- NS(id)
  plotlyOutput(ns("plot3"))
}

plot3_Server <- function(id, eventData2, rv, matrix, genesets){
  moduleServer(
    id,
    function(input, output, session){
      output$plot3 <- renderPlotly({
        req(eventData2())
        selected.gene.set <- rv$p2$x$data[[1]]$text[eventData2()]
        if(is(genesets(), "GeneSetCollection")){
          genes.toplot <- geneIds(genesets())[[selected.gene.set]]
        } else {
          genes.toplot <- genesets()[[selected.gene.set]]
        }
        mt <- match(genes.toplot, rownames(matrix()))
        x <-  matrix()[na.omit(mt), rv$sample.c]
        df <- as.data.frame(x)
        df$x <- as.numeric(df$x)
        df$Gene <- rownames(df)
        df$Sample <- rv$sample.c
        rv$p3 <- ggplot(data = df, aes(x=x, color = Sample, label = Gene)) +
          stat_density(geom="line", position = "identity") +
          geom_rug() + theme(legend.position = "none") +
          labs(x="Gene expression values in selected gene set and sample", y="Density") +
          xlim(as.numeric(range(matrix()))) +
          scale_color_manual("legend", values= rv$dd.col)
        ggplotly(rv$p3, tooltip = c("Gene", "x")) %>% style(hoverinfo="none", traces = 1) %>%
          layout(title = list(text = paste0('<br><sup><i>', selected.gene.set, '</i></sup>'),
                              font = list(size=15)))
      })
    }
  )
}
