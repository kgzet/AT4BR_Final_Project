# final project AT4BR
# Kinga Zajdel

# install.packages("reticulate")
# library(reticulate)
# source_python("FeLV_codon_usage_bias_counting.py")

library(shiny)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(plotly)

#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#

# data from .csv file from Python script 
table_of_codons <- read_csv("data/codon_usage.csv")
# sorting the table according to amino acid name
new_table_of_codons <- table_of_codons |> dplyr::arrange(aa, occurrence)
# preparing table with amino acids only, with no duplicates (for select input list)
choices <- distinct(new_table_of_codons, aa)

# Define server logic required to draw a bar plot
function(input, output, session) {
  
  # updateSelectInput to save the memory
  chosen_aa <- updateSelectInput(session ,"select_aa", "Chose an amino acid", choices$aa)
  
  # creating variables that will be used for plots
  x_axis_variable <- "nucleotide triplet(s) of"
  plot_title_numbers <- "bar plot of codon usage bias for triplet(s) coding"
  plot_title_percent <- "codon usage bias in percents for triplet(s) coding"
  
  
  output$distPlot <- renderPlotly({
    
    # which amino acid user wants to plot
    selected_aa <- input$select_aa
    # which kind of data user wants to plot (numbers of occurrence or percentage)
    name_of_data <- input$radio_data
    # table with only one amino acid, basis for the shown plot
    one_aa_table <- new_table_of_codons |> filter(aa == selected_aa)
    
    # a dictionary of amino acid shortcuts -> full names, for more clarity of the plots
    amino_name = switch(   
      selected_aa,   
      "A" = "Alanine",   
      "C" = "Cysteine",   
      "D" = "Aspartic acid",   
      "E" = "Glutamic acid", 
      "F" = "Phenylalanine", 
      "G" = "Glycine",
      "H" = "Histidine",
      "I" = "Isoleucine",
      "K" = "Lysine",
      "L" = "Leucine",
      "M" = "Methionine",
      "N" = "Asparagine",
      "P" = "Proline",
      "Q" = "Glutamine",
      "R" = "Arginine",
      "S" = "Serine",
      "T" = "Threonine",
      "V" = "Valine",
      "W" = "Tryptophan",
      "Y" = "Tyrosine",
      "*" = "stop"
    )
    
    # depending on user's choice one of the given plots will be shown
    # radio button ID: radio_data
    if(name_of_data == "numbers"){
      
      final_plot <- one_aa_table |>
        ggplot(one_aa_table, mapping = aes(x = codon, y = occurrence)) +    # using table with only one chosen aa
        geom_bar(stat = "identity", width = 0.8, col = input$radio_border, fill = input$radio_color) +  # user chooses decoration
        labs(x = paste(x_axis_variable, selected_aa), y = "occurrence", title = paste(plot_title_numbers, amino_name)) +
        # paste()allows to concatenate string variables, so labels will contain actual aa name 
        theme_minimal()
      
    }
    
    else if (name_of_data ==  "percentage"){
      
      # preparing data to count percentage
      # adding sum of all occurrences of chosen aa - mutate() allows to add a column to the tab
      percent_occurrence <- one_aa_table |>
        mutate(total_occurrence = sum(occurrence))
      
      # counting percentage of the chosen aa, and adding a column to the table
      percent_occurrence <- percent_occurrence |>
        mutate(percentage = (occurrence / total_occurrence) * 100)
      
      final_plot <- percent_occurrence |>     # using the table with only one chosen aa and added columns
        ggplot(percent_occurrence, mapping = aes(x = codon, y = percentage)) + 
        # geom_bar(stat = "identity", width = 0.8, col = "lightpink", fill = "orchid") +
        geom_bar(stat = "identity", width = 0.8, col = input$radio_border, fill = input$radio_color) +  # user chooses decoration
        labs(x = paste(x_axis_variable, selected_aa), y = "percentage (%)", title = paste(plot_title_percent, amino_name)) +
        # paste()allows to concatenate string variables, so labels will contain actual aa name 
        theme_minimal()
      
    }
    
    # printing one of the plots chosen by user - "if" statement decides which will be sent to the ggplotly() function
    ggplotly(final_plot)
    
  })
  
}

# No AI used.
