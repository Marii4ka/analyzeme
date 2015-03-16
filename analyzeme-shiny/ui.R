library(shiny)

shinyUI(
  fluidPage(
    titlePanel("analyzeMe"),
    mainPanel(
      p("Have you had your genome analyzed with 23andme? Have you wondered how this data can be analyzed? Well you're in the right place! Upload your raw data here, and we will analyze it and walk through what exactly it is we're looking at. "),
      p("Have you uploaded your file (or used the example)? Great! The analysis will take a while, but while you're waiting, we can talk about what we're doing!")
    )
  ) 
)