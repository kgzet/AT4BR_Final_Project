# final project AT4BR
# Kinga Zajdel

# install.packages("reticulate")
# library(reticulate)
# source_python("FeLV_codon_usage_bias_counting.py")

#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#

library(shiny)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(plotly)

table_of_codons <- read_csv("data/codon_usage.csv")
new_table_of_codons <- table_of_codons |> dplyr::arrange(aa, occurrence)
choices <- distinct(new_table_of_codons, aa)

# Define server logic required to draw a bar plot
function(input, output, session) {
      
      # selected_aa <- 'A'
      # col = input$radio_color
      # selected_aa <- updateSelectInput(session ,"select_aa", "Chose an amino acid", choices$aa)

      x_axis_variable <- "nucleotide triplets of "
      plot_title_numbers <- "bar plot of codon usage bias for amino acid:"
      plot_title_percent <- "codon usage bias in percents for amino acid "
      
      
      output$distPlot <- renderPlotly({
        
        selected_aa <- input$select_aa
        name_of_data <- input$radio_data
        # one_aa_table <- new_table_of_codons |> filter(aa == input$select_aa)
        one_aa_table <- new_table_of_codons |> filter(aa == selected_aa)
      
        
        if(name_of_data == "numbers"){
          
          final_plot <- one_aa_table |>
          # filter(aa==input$select_aa) |>
          ggplot(one_aa_table, mapping=aes(x = codon, y = occurrence)) +
          geom_bar(stat = "identity", width = 0.8, col = input$radio_border, fill = input$radio_color) +
          labs(x = paste(x_axis_variable, selected_aa), y = "occurrence", title = paste(plot_title_numbers, selected_aa)) +
          theme_minimal()
          
        }
        
        
        else if (name_of_data ==  "percentage"){
          
          percent_occurrence <- one_aa_table |>
            mutate(total_occurrence = sum(occurrence))
          
          percent_occurrence <- percent_occurrence |>
            mutate(percentage = (occurrence / total_occurrence) * 100)
          
          final_plot <- percent_occurrence |> 
            ggplot(percent_occurrence, mapping=aes(x = codon, y = percentage)) +
            # geom_bar(stat = "identity", width = 0.8, col = "lightpink", fill = "orchid") +
            geom_bar(stat = "identity", width = 0.8, col = input$radio_border, fill = input$radio_color) +
            labs(x = paste(x_axis_variable, selected_aa), y = "percentage (%)", title = paste(plot_title_percent, selected_aa)) +
            theme_minimal()
          
        }
        
        ggplotly(final_plot)

    })

}
