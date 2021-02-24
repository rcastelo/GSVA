gsva_information <- function(input, output, session, newY, genes) {
  
  
  #Rendering text1
  output$text1 <- renderUI({
    HTML(paste("<br/>", "\t To see the Empirical Cumulative Distribution Function 
    of a Sample, click on its line in this plot and go
      to the 'Gene.Set' Panel", "<br/>", sep="<br/>"))
  })
  
  # Rendering graph1
  dat.t <- melt(as.data.table(generated_gsva, keep.rownames = "gene.sets"), 
                variable.name = "Sample", id.vars="gene.sets")
  n <- length(levels(dat.t$Sample))
  dd.col <- hcl(h = seq(15, 375, length=n), l = 65, c = 100)[1:n]
  names(dd.col)  <- levels(dat.t$Sample)
  
  
  output$plot <- renderPlotly({
    p <- ggplot(data = dat.t, aes(x=value, color=Sample)) +
      stat_density(geom="line", position = "identity") +
      theme(legend.position = "none") + labs(x="GSVA Scores", y="Density") +
      scale_color_manual("Legend", values = dd.col)
    ggplotly(p, tooltip = "Sample", source = "click1")
  }) 
  
  # Rendering table
  resultInformation <- matrix(data = c(nrow(generated_gsva),
                                       ncol(generated_gsva)),
                              nrow = 1, ncol = 2)
  colnames(resultInformation) <- c("Nr. of gene sets", "Nr. of samples")
  output$result <- renderTable(resultInformation)
  
  #Rendering text2
  output$text2 <- renderUI({
    title1 <- sample.c()
    h2(tags$b(title1), align ="center")
  })
  
  #Rendering text3
  output$text3 <- renderUI({
    HTML(paste("<br/>", "\t To see the Kernel Density Estimation of genes of 
    any given Gene Set in this Sample,  click on any point in this plot and a
    second plot will appear bellow it", "<br/>", sep="<br/>"))
  })
  
  #Rendering graph2
  eventData1 <- reactive({
    event_data("plotly_click", source = "click1")
  })
  
  sample.c <- reactive({
    req(eventData1())
    ind <- eventData1()$curveNumber+1
    colnames(generated_gsva)[ind]
  })
  
  plot2 <- reactive({
    req(sample.c())
    data <- dat.t[Sample==sample.c()]
    p <- ggplot(data = data, aes(x=value, color=Sample)) +
      stat_ecdf(geom="point") + theme(legend.position = "none") + 
      labs(x="GSVA Scores in selected sample", y="Empirical Cumulative Density") +
      scale_color_manual("Legend", values = dd.col)
    p <- ggplotly(p, source="click2") %>% style(text=data$gene.sets)
  })
  
  output$plot2 <- renderPlotly({
    req(plot2())
    plot2()
  })
  
  # Rendering graph 3
  eventData2 <- reactive({
    event_data("plotly_click", source = "click2")
  })
  
  gene.set <- reactive({
    plot2()$x$data[[1]]$text[eventData2()$pointNumber+1]
  })
  
  output$plot3 <- renderPlotly({
    req(eventData2())
    if(is(genes, "GeneSetCollection")){
      genes.toplot <- geneIds(genes)[[gene.set()]]  
    } else {
      genes.toplot <- genes[[gene.set()]]
    }
    mt <- match(genes.toplot, rownames(newY))
    x <-  newY[na.omit(mt), sample.c()]
    df <- as.data.frame(x)
    df$x <- as.numeric(df$x)
    df$Gene <- rownames(df)
    df$Sample <- sample.c()
    p1 <- ggplot(data = df, aes(x=x, color = Sample, label = Gene)) +
      stat_density(geom="line", position = "identity") +
      geom_rug() + theme(legend.position = "none") + 
      labs(x="Gene Expressions in selected sample", y="Density") +
      xlim(as.numeric(range(newY))) +
      scale_color_manual("legend", values= dd.col)
    ggplotly(p1, tooltip = c("Gene", "x")) %>% style(hoverinfo="none", traces = 1) %>%
      layout(title = list(text = paste0('<br><sup><i>', gene.set(), '</i></sup>'),
                          font = list(size=15)))
  })
  
  # Rendering Session Info
  output$sessionInfo <- renderPrint({
    sessionInfo()
  })
  
  download_handler(input,output,session)
  
}