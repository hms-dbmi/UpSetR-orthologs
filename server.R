library(shiny)
library(UpSetR)
library(biomaRt)
library(data.table)
source("ensembl.R")
shinyServer(function(input, output, session){
  options(shiny.maxRequestSize=30*1024^2)
  datasets <- reactive({
    names <- read.csv("genes.csv")
    names <- sort(names$description)
    return(names)
  })
  
  
  output$genes <- renderUI({
    genes <- selectInput('Select', h6("Select organisms:"),
                         choices = as.character(datasets()),
                         multiple = T, selectize = T, selected = NULL)
    return(genes)
  })
  
    organisms <- reactive({
      organisms <- as.character(c(input$Select))
    })
  
  organisms <- eventReactive(input$goButton,
{
  as.character(input$Select)
})

  species <- reactive({
    species <- read.csv("genes.csv")
    species <- species[which(species$description %in% organisms()),]
    species <- head(as.character(species$dataset), 5)
  })
  
  data <- reactive({
   m <- (length(species())/10)
      if(length(species()) > 6){warning("6 species maximum")}
      else{
        datasets <- list()
        withProgress(message = 'Selecting data sets from Biomart', value = 0, max = m,
                     detail = "Data set 0", {
          Sys.sleep(0.1)
        for(i in seq_along(species())){
          incProgress(0.1, detail = paste("Data set", i))
          datasets[[i]] <- getEnsemblMarts(species()[i])
        }
       })
      }
      return(datasets)
  })
  
  namesData <- reactive({
   namesData <- getNames(species())
   return(namesData)
  })
  
  homologs <- reactive({
    homologs <- list()
    m <- (length(species())/10)
    withProgress(message = 'Querying data sets from BioMart Website', value = 0, max =m,
                 detail= 'Querying on data set 0', {
                   Sys.sleep(0.1)
      for(i in seq_along(data())){
        incProgress(0.1, detail = paste("Querying on data set", i))
        homologs[[i]] <- getBM(attributes = c("ensembl_gene_id",namesData()[[i]]), mart = data()[[i]])
        new_colname <- paste(unlist(strsplit(species()[i], split = "_"))[1], "_homolog_ensembl_gene", sep = "", collapse="")
        colnames(homologs[[i]])[1] <- new_colname
        homologs[[i]] <- homologs[[i]][ ,order(colnames(homologs[[i]]))]
      }
    })
      alldata <- rbindlist(homologs)
      alldata <- data.frame(alldata[!duplicated(alldata),])
      return(alldata)
  })
  
  upsetdata <- reactive({
    withProgress(message = 'Converting to upset data', value = 0, {
      upsetdata <- BioMartConverter(homologs())
      setProgress(1)
      return(upsetdata)
    })
  })
  
  
  My_data <- reactive({  
    My_data <- list()
    inFile <- input$files
    if (is.null(inFile))
      return(NULL)
    else{
    My_data <- lapply(inFile$datapath, read.csv)
    }
return(My_data)
  })
  
  upsetdata_manual <- reactive({
    upsetdata_manual <- orthologs_manual(My_data())
    upsetdata_manual <- ManualBioMartConverter(upsetdata_manual)
    return(upsetdata_manual)
  })


output$plot_text <- renderUI({
  if(length(input$Select) == 0){
  text1 <- "This is where your plot will be displayed."
  text2 <- "Some datasets can be extremely large, so plotting may take some time. Please be patient."
  text3 <- "It is best to refresh after each run of the plot."
  HTML(paste(text1, text2, text3, sep = '<br/><br/>'))
  }
  else{
    return()
  }
})

output$mplot_text <- renderUI({
  if(is.null(My_data())){
    text <- "This is where your plot will be displayed!"
    HTML(text)
  }
  else{
    return()
  }
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
    upset(upsetdata(), order.matrix="freq", empty.intersections = "on")
  
  dev.off()
  # Return a list
  list(src = outfile,
       width = width,
       height = height)
}, deleteFile = TRUE)

output$down <- downloadHandler(
  
  filename = function(){
    paste("UpSetR-orthologs", input$filetype, sep =".")
  }, 
  content = function(file){
    width  <- session$clientData$output_plot_width
    height <- ((session$clientData$output_plot_height))
    pixelratio <- session$clientData$pixelratio
    if(input$filetype == "png")
      png(file, width=width*pixelratio, height=height*pixelratio,
          res=72*pixelratio)
    else
      pdf(file,width = 22, height = 14)
    upset(upsetdata(), order.matrix="freq", empty.intersections = "on")
    
    dev.off()
  }
)

output$manplot <- renderImage({

  width  <- session$clientData$output_manplot_width
  height <- ((session$clientData$output_manplot_height) * 1.7)
  pixelratio <- session$clientData$pixelratio
  # A temp file to save the output. It will be deleted after renderImage
  # sends it, because deleteFile=TRUE.
  
  if(length(My_data()) == 0){
    stop()
  }
  
  outfile <- tempfile(fileext='.png')
  
  # Generate a png
  if(!is.null(My_data())){
  png(outfile, width=width*pixelratio, height=height*pixelratio,
      res=72*pixelratio)
    upset(upsetdata_manual(), order.matrix="freq", empty.intersections = "on")
  }
  dev.off()
  # Return a list
  list(src = outfile,
       width = width,
       height = height)
}, deleteFile = TRUE)


output$mandown <- downloadHandler(
  
  filename = function(){
    paste("UpSetR-file", input$ftype, sep =".")
  }, 
  content = function(file){
    width  <- session$clientData$output_manplot_width
    height <- ((session$clientData$output_manplot_height))
    pixelratio <- session$clientData$pixelratio
    if(input$ftype == "png")
      png(file, width=width*pixelratio, height=height*pixelratio,
          res=72*pixelratio)
    else
      pdf(file,width = 22, height = 14)
    upset(upsetdata_manual(), order.matrix="freq", empty.intersections = "on")
    
    dev.off()
  }
)

})