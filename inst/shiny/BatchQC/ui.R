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
shinyUI(navbarPage("BatchQC", id="BatchQC", fluid=TRUE, 
                   tabPanel("Differential Expression",
                            sidebarLayout(
                              sidebarPanel(
                                numericInput('ncSamples', 'No. of Sample(s) Per Condition', 
                                             if (maxcondElems>defaultDisp) defaultDisp else maxcondElems,
                                             min = 1, max = maxcondElems),
                                numericInput('noSamples', 'No. of Sample(s) Per Batch', 
                                             if (maxbatchElems>defaultDisp) defaultDisp else maxbatchElems,
                                             min = 1, max = maxbatchElems),
                                checkboxInput("sortbybatch", "Sort By Batches First", FALSE),
                                checkboxInput("colbybatch", "Color By Batches", FALSE),
                                checkboxInput("combatDE", "ComBat", FALSE),
                                # This makes web page load the JS file in the HTML head.
                                #singleton(
                                #  tags$head(tags$script(src = "message-handler.js"))
                                #)
                                singleton(
                                  tags$head(tags$script(
                                    'Shiny.addCustomMessageHandler("testmessage",
                                    function(message) {
                                        alert(message);
                                      }
                                    );'
                                  ))
                                )
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
                              checkboxInput("combatHM", "ComBat", FALSE)
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
                                checkboxInput("combatCD", "ComBat", FALSE)
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
                          checkboxInput("combatPCA", "ComBat", FALSE)
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
                                             min = 1, max = nrow(shinyInput$delta.hat)),
                                br(),
                                actionButton("runCombat", "Run ComBat"),
                                p("Click the button to run ComBat and see the results in the other tabs.")
                              ),
                              mainPanel(
                                tabsetPanel(
                                  id="CombatMain",
                                  tabPanel("ComBat Plots", plotOutput("densityQQPlots")),
                                  tabPanel("Summary", verbatimTextOutput("kstest")),
                                  tabPanel("ComBat Output", value="ComBatOutput", verbatimTextOutput("combatOutText"))
                                )
                              )
                            )
                   )
      
))