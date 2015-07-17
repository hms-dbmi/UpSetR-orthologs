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
                                width = 3
                              ),
                              mainPanel(
                                imageOutput("plot"), width = 9
                              )
                            )
                   )
                   
))