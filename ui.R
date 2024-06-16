# final project AT4BR
# Kinga Zajdel

# install.packages("reticulate")
# library(reticulate)

#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#

library(shiny)
library(plotly)

# Define UI for application that draws a bar plot
fluidPage(

    # Application title
    titlePanel("Codon bias on the example of FeLV"),

    # Sidebar with a slider input for number of bins
    fluidRow(
      
      sidebarPanel(
        selectInput("select_aa", "Select an amino acid:", choices=""),
        # selectInput(inputId = "select_aa", label = "Select an amino acid:", choices = choices),
        radioButtons(inputId = "radio_color", label = "Choose color:", choices = c("grey", "darkred", "violet", "lightblue")),
        radioButtons(inputId = "radio_border", label = "Choose border:", choices = c("black", "sienna", "plum4", "lightpink")),
        radioButtons(inputId = "radio_data", label = "Choose data to show:", choices = c("numbers", "percentage")),
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        plotlyOutput("distPlot")
      )
      
    )
)

# No AI used.
