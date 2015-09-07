library(shiny)
library(ggvis)
library(d3heatmap)

shinyUI(navbarPage("BatchQC",
                   tabPanel("Box Plots",
                            sidebarLayout(
                              sidebarPanel(
                                numericInput('noSamples', 'Number of Samples', 1,
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
                   tabPanel("PCA Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          numericInput('xcol', 'Principal Component (x-axis)', 1,
                            min = 1, max = 50),
                          numericInput('ycol', 'Principal Component (y-axis)', 2,
                            min = 1, max = 50)
                        ),
                        mainPanel(
                          tabsetPanel(
                          tabPanel("PCA", ggvisOutput("plot")),
                            tabPanel("SVD", plotOutput("svd")),
                            tabPanel("Summary", verbatimTextOutput("PCAsummary")),
                            tabPanel("Table",tableOutput("PCAtable"))
                          )
                        )
                      )
                    ),
                   tabPanel("Outliers", plotOutput("outliers")),
                   tabPanel("Heatmaps",
                     tabsetPanel(
                       tabPanel("Heatmap",
                                #selectInput("palette", "Palette", c("RdBu","Greens", "Blues")),
                                checkboxInput("cluster1", "Apply clustering"),
                                d3heatmapOutput("heatmap")
                        ),
                       tabPanel("Sample Correlations",
                               checkboxInput("cluster2", "Apply clustering"),
                               d3heatmapOutput("correlation")
                       )
                    )
                  ),
                   tabPanel("Combat",
                            sidebarLayout(
                              sidebarPanel(
                                numericInput('batches', 'Batch', 1,
                                             min = 1, max = nrow(shinyInput$delta.hat))
                              ),
                              mainPanel(
                                tabsetPanel(
                                  tabPanel("Combat Plots", plotOutput("densityQQPlots")),
                                  tabPanel("Summary", verbatimTextOutput("kstest"))
                                )
                              )
                            )
                   )
))