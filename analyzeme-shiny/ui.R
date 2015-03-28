library(shiny)

shinyUI(
  fluidPage(
    fluidRow(
      column(
        12,
        offset=1,
        h1("analyzeMe"),
        h4("by Arun Durvasula")
      )
    ),
    fluidRow(
      column(
          2,
          offset=1,
          h3("About"),
          p("Built with shiny and R")
        ),
      column(
        8,
        p("Have you had your genome analyzed with 23andme? Have you wondered how this data can be analyzed? Well you're in the right place! Upload your raw data here, and we will analyze it and walk through what exactly it is we're looking at. "),
        p("Privacy: your data will be deleted in a week and will never be looked at."),
        fileInput(inputId="inFile", label="Please select a zipped 23andme raw data file.", accept=c("application/x-gzip")),
        p("or"),
        checkboxInput(inputId="useExample", label="Use Example", value=T),
        p("Have you uploaded your file (or used the example)? Great! The analysis will take a while, but while you're waiting, let's talk about what we're doing."),
        h2("The data"),
        p("The data you receive from 23 and Me is the result of a genotyping experiment. This is not the same as sequencing, but when you are only trying to find out what mutations exist at single base pairs, it is a much cheaper alternative to sequencing."),
        p("What's a base pair? Well let's start at the beginning. Our bodies are made up of cells. These cells come in different types (for example: skin cells, brian cells, liver cells) and each one contains several molecules of DNA that make up your genome. Your genome has all the information needed for your body to create proteins (and do other things too). Proteins carry out important functions in your cells like transporting nutrients in and out, breaking down compounds for energy, and replicating your DNA to make more cells."), 
        p("Sometimes diseases are caused by mutations in genes that create proteins."),

        plotOutput("snpsByChr"),
        dataTableOutput('riskTable'),
        conditionalPanel(
          condition = "output.useExample == true",
          imageOutput("allPCA.image")
          ),
        conditionalPanel(
          condition = "output.useExample == false",
          plotOutput("allPCA.plot")
        )
      )
    )
  ) 
)