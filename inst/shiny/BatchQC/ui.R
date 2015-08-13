library(shiny)
library(ggvis)
library(d3heatmap)

shinyUI(navbarPage("BatchQC",
                   tabPanel("Box Plots",
                            sidebarLayout(
                              sidebarPanel(
                                numericInput('noGenes', 'Number of Genes', 10,
                                             min = 1, max = ncol(data.matrix))
                              ),
                              mainPanel(
                                tabsetPanel(
                                  tabPanel("Boxplot",ggvisOutput("Boxplot")), 
                                  tabPanel("Summary", verbatimTextOutput("BPsummary")),
                                  tabPanel("Table", tableOutput("BPtable")) 
                                )
                              )
                            )
                   ),
                   tabPanel("PCA",
                      sidebarLayout(
                        sidebarPanel(
                          numericInput('xcol', 'Principal Component (x-axis)', 1,
                            min = 1, max = 3),
                          numericInput('ycol', 'Principal Component (y-axis)', 2,
                            min = 1, max = 3)
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Plot", ggvisOutput("plot")), 
                            tabPanel("Summary", verbatimTextOutput("PCAsummary")),
                            tabPanel("Table",tableOutput("PCAtable"))
                          )
                        )
                      )
                    ),
                   tabPanel("Heatmap",
                            selectInput("palette", "Palette", c("RdBu","Greens", "Blues")),
                            checkboxInput("cluster", "Apply clustering"),
                            d3heatmapOutput("heatmap")
                   ),
                   tabPanel("Combat",
                            sidebarLayout(
                              
                              sidebarPanel(
                                numericInput('batches', 'Batch', 1,
                                             min = 1, max = 2)#nrow(gamma.hat))
                              ),
                              
                              mainPanel(
                                tabsetPanel(
                                  tabPanel("Plot", ggvisOutput("density")), 
                                  tabPanel("Summary")
                                )
                              )
                            )
                   )
              
                          
))

