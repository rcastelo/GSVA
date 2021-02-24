fluidPage( 
  use_waiter(),
  use_hostess(),
  theme = shinytheme("spacelab"),	
  tags$head(tags$style("#errorsGsva{color: red;
                                 font-size: 20px;
                                 font-style: italic;
                                 }"
                       )),
  titlePanel(
    fluidRow(
      column(6,
        h2("GSVA Shiny App", align="left")),
      column(6,
        tags$img(src="GSVA.png", align="right", height=75, width=75))
      ), windowTitle="GSVA"),
	fluidRow(
	  selectDataInput("dataInput"),
	  mainDataInput("mainInput"),
	  argumentsDataInput("argumentsInput")
	)
)
