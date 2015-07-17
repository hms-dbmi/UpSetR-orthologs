library(shiny)
library(UpSetR)
library(biomaRt)
source("ensembl.R")
shinyServer(function(input, output, session){
  
  datasets <- reactive({
    names <- listDatasets(useMart("ensembl"))
    names <- sort(names$description)
    return(names)
  })
  
  output$genes <- renderUI({
    genes <- selectInput('Select', h6("Select organisms:"),
                         choices = as.character(datasets()),
                         multiple = T, selectize = T, selected = NULL)
    return(genes)
  })
  
  #   organisms <- reactive({
  #     organisms <- as.character(c(input$Select))
  #   })
  
  organisms <- eventReactive(input$goButton,
{
  as.character(input$Select)
})

output$plot <- renderImage({
  
  width  <- session$clientData$output_plot_width
  height <- ((session$clientData$output_plot_height)*1.7)
  pixelratio <- session$clientData$pixelratio
  # A temp file to save the output. It will be deleted after renderImage
  # sends it, because deleteFile=TRUE.
  outfile <- tempfile(fileext='.png')
  
  # Generate a png
  png(outfile, width=width*pixelratio, height=height*pixelratio,
      res=72*pixelratio)
  if(!is.null(organisms()))
    UpSetRensembl(organisms())
  
  dev.off()
  # Return a list
  list(src = outfile,
       width = width,
       height = height)
}, deleteFile = TRUE)

})