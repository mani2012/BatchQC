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
nbatch <- nlevels(as.factor(shinyInput$batch))
shinyUI(navbarPage("BatchQC", id="BatchQC", fluid=TRUE, 
                   tabPanel("Summary",
                            tabsetPanel(
                              tabPanel("Confounding", 
                                       br(),
                                       h3("Number of samples in each Batch and Condition"),
                                       tableOutput("CondBatchTable"),
                                       br(),
                                       h3("Measures of confounding between Batch and Condition"),
                                       tableOutput("ConfoundTable")
                              ),
                              tabPanel("Variation Analysis", 
                                       radioButtons('batchVA', 'Batch Adjustment',
                                                    c('None'=0, 'Combat'=1,'SVA'=2), 0),
                                       br(),
                                       h3("Variation explained by Batch and Condition"),
                                       plotOutput("VariationPlot"),
                                       br(),
                                       tableOutput("VariationTable")
                              ),
                              tabPanel("P-value Analysis", 
                                       radioButtons('batchPA', 'Batch Adjustment',
                                                    c('None'=0, 'Combat'=1,'SVA'=2), 0),
                                       br(),
                                       h3("Distribution of Batch and Condition Effect p-values Across Genes"),
                                       tableOutput("PvalueTable"),
                                       br(),
                                       plotOutput("BatchPvaluePlot"),
                                       br(),
                                       plotOutput("ConditionPvaluePlot")
                              )
                            )
                   ),
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
                                radioButtons('batchDE', 'Batch Adjustment',
                                             c('None'=0, 'Combat'=1,'SVA'=2), 0),
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
                                  tabPanel("Expression Plots",ggvisOutput("DiffExPlot")), 
                                  tabPanel("Summary", verbatimTextOutput("DEsummary")),
                                  tabPanel("Table", tableOutput("DEtable")), 
                                  tabPanel("LIMMA",tableOutput("LimmaTable")),
                                  tabPanel("GLS",tableOutput("GlsTable")), 
                                  tabPanel("Mixed Effects",tableOutput("MixEffTable")) 
                                )
                              )
                            )
                   ),
                   tabPanel("Median Correlations", plotOutput("outliers")),
                   tabPanel("Heatmaps",
                            tabsetPanel(
                              tabPanel("Heatmap",
                                       #selectInput("palette", "Palette", c("RdBu","Greens", "Blues")),
                                       radioButtons('batchHM1', 'Batch Adjustment',
                                                    c('None'=0, 'Combat'=1,'SVA'=2), 0),
                                       checkboxInput("cluster1", "Apply clustering"),
                                       d3heatmapOutput("heatmap")
                              ),
                              tabPanel("Sample Correlations",
                                       radioButtons('batchHM2', 'Batch Adjustment',
                                                    c('None'=0, 'Combat'=1,'SVA'=2), 0),
                                       checkboxInput("cluster2", "Apply clustering"),
                                       d3heatmapOutput("correlation")
                              )
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
                                checkboxInput("colbybatchCD", "Color By Batches", FALSE),
                                radioButtons('batchCD', 'Batch Adjustment',
                                             c('None'=0, 'Combat'=1,'SVA'=2), 0)
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
                          checkboxInput("colbybatchPCA", "Color By Batches", FALSE),
                          radioButtons('batchPCA', 'Batch Adjustment',
                                       c('None'=0, 'Combat'=1,'SVA'=2), 0)
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("PCA", ggvisOutput("plot")),
                            tabPanel("Summary", verbatimTextOutput("PCAsummary")),
                            tabPanel("Table",tableOutput("PCAtable")),
                            tabPanel("Explained Variation",tableOutput("PCAExplainedVariation"))
                          )
                        )
                      )
                    ),
                   tabPanel("Shape",
                            sidebarLayout(
                              sidebarPanel(
                                numericInput('batchnum', 'Batch', 1,
                                             min = 1, max = nbatch),
                                radioButtons('batchShape', 'Batch Adjustment',
                                             c('None'=0, 'Combat'=1,'SVA'=2), 0)
                              ),
                              mainPanel(
                                tabsetPanel(
                                  tabPanel("Manova", plotOutput("Manova")),
                                  tabPanel("SO Plot", plotOutput("SOplot"))
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
                                radioButtons('optionMeanOnly', 'Batch Adjustment Option',
                                             c('Mean and Variance'=0, 'Mean only'=1), 0),
                                br(),
                                radioButtons('optionParameter', 'Batch Parameter Distribution',
                                             c('Parametric'=0, 'Non-Parametric'=1), 0),
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
                   ),
                   tabPanel("SVA",
                            sidebarLayout(
                              sidebarPanel(
                                checkboxInput("fsvaOption", "Frozen SVA (Default: Regression Adjusted)", FALSE),
                                actionButton("runSVA", "Run SVA"),
                                p("Click the button to run SVA and see the results in the other tabs.")
                              ),
                              mainPanel(
                                tabsetPanel(
                                  id="SVAMain",
                                  tabPanel("Summary", verbatimTextOutput("svasummary")),
                                  tabPanel("SVA Output", value="SVAOutput", verbatimTextOutput("SVAOutText"))
                                )
                              )
                            )
                   )
                   
))