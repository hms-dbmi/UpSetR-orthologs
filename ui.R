library(shiny)

shinyUI(
  navbarPage("UpSetR - Ensembl",
             tabPanel("Instructions",
                      mainPanel()     
             ),
             tabPanel("Ensembl via R",
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
             ),
             tabPanel("Ensembl via files",
                      sidebarLayout(
                        sidebarPanel(
                        fluidRow(
                          fileInput('files', 'Upload files', multiple = TRUE, accept = c(
                            'text/csv', 'text/comma-separated-values', 'text/tab-separated-values', '.csv', '.tsv')),
                          br(),
                          radioButtons(inputId = "ftype", label = "File type", choices = list("png", "pdf")),
                          br(),
                          downloadButton(outputId = "mandown", label = "Download!")
                        )
                        ,width =3),
                      mainPanel(
                        htmlOutput("mplot_text"),
                        imageOutput('manplot'),
                      width = 9)
             )
             )
             
  ))