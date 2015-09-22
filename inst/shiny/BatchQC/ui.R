library(shiny)
library(ggvis)
library(d3heatmap)

minbatch <- function(batch1){
  batch2 <- as.factor(batch1)
  batch3 <- split(batch1,batch2)
  return(min(unlist(lapply(1:length(batch3), function(x) length(batch3[[x]])))))
}

maxbatchElems <- minbatch(shinyInput$batch)
maxcondElems <- minbatch(shinyInput$condition)
defaultDisp <- 30
shinyUI(navbarPage("BatchQC",
                   tabPanel("Box Plots",
                            sidebarLayout(
                              sidebarPanel(
                                numericInput('noSamples', 'No. of Sample(s) Per Batch', 
                                             if (maxbatchElems>defaultDisp) defaultDisp else maxbatchElems,
                                             min = 1, max = maxbatchElems),
                                checkboxInput("combat", "ComBat", FALSE)
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
                   tabPanel("Differential Expression",
                            sidebarLayout(
                              sidebarPanel(
                                numericInput('ncSamples', 'No. of Sample(s) Per Condition', 
                                             if (maxcondElems>defaultDisp) defaultDisp else maxcondElems,
                                             min = 1, max = maxcondElems),
                                checkboxInput("colbybatch", "Color By Batches", FALSE),
                                checkboxInput("combat", "ComBat", FALSE)
                              ),
                              mainPanel(
                                tabsetPanel(
                                  tabPanel("Differential Expression",ggvisOutput("DiffExPlot")), 
                                  tabPanel("Summary", verbatimTextOutput("DEsummary")),
                                  tabPanel("Table", tableOutput("DEtable")) 
                                )
                              )
                            )
                   ),
                   tabPanel("Median Correlations", plotOutput("outliers")),
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
                              ),
                              checkboxInput("combat", "ComBat", FALSE)
                            )
                   ),
                   tabPanel("Circular Dendrogram",
                            sidebarLayout(
                              sidebarPanel(
                                selectInput("CorMethod", "Correlation Method:",
                                            c("Spearman"="spearman")),
                                selectInput("AggMethod", "Agglomeration Method:",
                                            c("Complete"="complete",
                                              "Ward" = "ward.D2",
                                              "Average" = "average",
                                              "McQuitty" = "mcquitty",
                                              "Single" = "single")),
                                checkboxInput("combat", "ComBat", FALSE)
                              ),
                              mainPanel(
                                plotOutput("circos", width = "100%")
                              )
                            )
                   ),
                   tabPanel("PCA Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          numericInput('xcol', 'Principal Component (x-axis)', 1,
                            min = 1, max = 50),
                          numericInput('ycol', 'Principal Component (y-axis)', 2,
                            min = 1, max = 50),
                          checkboxInput("combat", "ComBat", FALSE)
                        ),
                        mainPanel(
                          tabsetPanel(
                          tabPanel("PCA", ggvisOutput("plot")),
                            tabPanel("Summary", verbatimTextOutput("PCAsummary")),
                            tabPanel("Table",tableOutput("PCAtable"))
                          )
                        )
                      )
                    ),
                   tabPanel("ComBat",
                            sidebarLayout(
                              sidebarPanel(
                                numericInput('batches', 'Batch', 1,
                                             min = 1, max = nrow(shinyInput$delta.hat))
                              ),
                              mainPanel(
                                tabsetPanel(
                                  tabPanel("ComBat Plots", plotOutput("densityQQPlots")),
                                  tabPanel("Summary", verbatimTextOutput("kstest"))
                                )
                              )
                            )
                   )
      
))