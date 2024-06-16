# final project AT4BR
# Kinga Zajdel

# install.packages("reticulate")
# library(reticulate)

library(shiny)
library(plotly)

#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#

# Define UI for application that draws a bar plot
fluidPage(

    # Application title
    titlePanel("Codon bias on the example of FeLV"),

    # Sidebar with a slider input for number of bins
    fluidRow(
      
      sidebarPanel(
        
        selectInput("select_aa", "Select an amino acid:", choices=""),  # updateSelectInput is more memory-saving
        # selectInput(inputId = "select_aa", label = "Select an amino acid:", choices = choices), # standard selectInput
        
        # the radio button to choose kind of data: number of occurrences or percentage
        radioButtons(inputId = "radio_data", label = "Choose data to show:", choices = c("numbers", "percentage")),
        
        # some radio buttons to choose the decoration of the plot
        radioButtons(inputId = "radio_color", label = "Choose color:", choices = c("grey", "darkred", "violet", "lightblue")),
        radioButtons(inputId = "radio_border", label = "Choose border:", choices = c("black", "sienna", "plum4", "lightpink")),
        
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        plotlyOutput("distPlot")
      )
      
    )
)

# No AI used.
