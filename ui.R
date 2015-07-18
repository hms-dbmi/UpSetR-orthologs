library(shiny)

shinyUI(navbarPage("UpSetR - Ensembl",
                   tabPanel("Instructions",
                            mainPanel()     
                   ),
                   tabPanel("UpSetR Plot",
                            sidebarLayout(
                              sidebarPanel(
                                htmlOutput("genes"),
                                actionButton("goButton", "Submit Genes"),
                                br(),br(),
                                radioButtons(inputId = "filetype", label = "File type", choices = list("png", "pdf")),
                                downloadButton(outputId = "down", label = "Download!"),
                                br(), br(),
                                tags$a(href = "https://upsetr.shinyapps.io/UpSetR-orthologs", "Refresh"),
                                width = 3
                              ),
                              mainPanel(
                                htmlOutput("plot_text"),
                                imageOutput("plot"), width = 9
                              )
                            )
                   )
                   
))