library(shiny)
library(pander)
library(ggvis)
library(d3heatmap)

minbatch <- function(batch1){
    batch2 <- as.factor(batch1)
    batch3 <- split(batch1,batch2)
    return(min(unlist(lapply(1:length(batch3), 
        function(x) length(batch3[[x]])))))
}

shinyInput <- getShinyInput()
maxbatchElems <- minbatch(shinyInput$batch)
maxcondElems <- minbatch(shinyInput$condition)
defaultDisp <- 30
defaultTopGenesDisp <- 10
defaultGenesDisp <- 50
maxGenes <- dim(shinyInput$data)[1]
nbatch <- nlevels(as.factor(shinyInput$batch))
shinyUI(navbarPage("BatchQC", id="BatchQC", fluid=TRUE, 
    tabPanel("Summary",
        sidebarLayout(
            sidebarPanel(
                radioButtons('batchVA', 'Batch Adjustment',
                    c('None'=0, 'Combat'=1,'SVA'=2), 0),
                sliderInput('noGenesVA', 
                    'No. of Genes to include in the table', value =
                    if (maxGenes>defaultGenesDisp) defaultGenesDisp 
                    else maxGenes, min = 1, max = maxGenes, step=1),
                width=3
            ),
            mainPanel(
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
                        h3("Variation explained by Batch and Condition"),
                        plotOutput("VariationPlot"),
                        br(),
                        tableOutput("VariationTable")
                    ),
                    tabPanel("P-value Analysis", 
                        h3(
            "Distribution of Batch and Condition Effect p-values Across Genes"),
                        tableOutput("PvalueTable"),
                        plotOutput("BatchPvaluePlot"),
                        br(),
                        plotOutput("ConditionPvaluePlot")
                    )
                ), width=9
            )
        )
    ),
    tabPanel("Differential Expression",
        sidebarLayout(
            sidebarPanel(
                sliderInput('ncSamples', 'No. of Sample(s) Per Condition', 
                    value=if (maxcondElems>defaultDisp) defaultDisp 
                    else maxcondElems, min = 1, max = maxcondElems,step=1),
                sliderInput('noSamples', 'No. of Sample(s) Per Batch', 
                    value=if (maxbatchElems>defaultDisp) defaultDisp 
                    else maxbatchElems, min = 1, max = maxbatchElems,step=1),
                checkboxInput("sortbybatch", 
                    "Sort By Batch First (Default: Sort By Condition First)", 
                    FALSE),
                checkboxInput("colbybatch", 
                    "Color By Batch (Default: Color By Condition)", FALSE),
                radioButtons('batchDE', 'Batch Adjustment',
                    c('None'=0, 'Combat'=1,'SVA'=2), 0),
                sliderInput('noGenes', 
                    'No. of top Differentially Expressed Genes to display', 
                    value=if (maxGenes>defaultTopGenesDisp) defaultTopGenesDisp 
                    else maxGenes, min = 1, max = maxGenes, step=1),
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
                ),
                width=3
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Expression Plots",ggvisOutput("DiffExPlot")), 
                    tabPanel("Summary", verbatimTextOutput("DEsummary")),
                    #tabPanel("Table", tableOutput("DEtable")), 
                    tabPanel("LIMMA",tableOutput("LimmaTable"))
                    #tabPanel("GLS",tableOutput("GlsTable")), 
                    #tabPanel("Mixed Effects",tableOutput("MixEffTable"))
                ), width=9
            )
        )
    ),
    tabPanel("Median Correlations", 
        sidebarLayout(
            sidebarPanel(
                radioButtons('batchMC', 'Batch Adjustment',
                    c('None'=0, 'Combat'=1,'SVA'=2), 0),
                width=3
            ),
            mainPanel(
                plotOutput("outliers", width="100%"), width=9)
        )
    ),
    tabPanel("Heatmaps",
        sidebarLayout(
            sidebarPanel(
                #selectInput("palette", "Palette", c("RdBu","Greens", "Blues")),
                radioButtons('batchHM', 'Batch Adjustment',
                    c('None'=0, 'Combat'=1,'SVA'=2), 0),
                checkboxInput("cluster", "Apply clustering"),
                sliderInput('noGenesHM', 
                    'No. of Genes to include in heatmap display', value =
                    if (maxGenes>defaultGenesDisp) defaultGenesDisp 
                    else maxGenes, min = 1, max = maxGenes, step=1),
                width=3
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Heatmap",
                        d3heatmapOutput("heatmap", width="100%", height="550px")
                    ),
                    tabPanel("Sample Correlations",
                        d3heatmapOutput("correlation", width="100%", 
                            height="550px")
                    )
                ) , width=9
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
                checkboxInput("colbybatchCD", 
                    "Color By Batch (Default: Color By Condition)", FALSE),
                radioButtons('batchCD', 'Batch Adjustment',
                    c('None'=0, 'Combat'=1,'SVA'=2), 0),
                width=3
            ),
            mainPanel(
                plotOutput("circos", height = "550px", width = "100%"),
                width=9
            )
        )
    ),
    tabPanel("PCA Analysis",
        sidebarLayout(
            sidebarPanel(
                selectInput('xcoln', 'Principal Component (x-axis)', choices=
                    as.character(1:(if (is.null(shinyInput$pc)) 1 
                    else ncol(shinyInput$pc))), selected="1"),
                selectInput('ycoln', 'Principal Component (y-axis)', choices=
                    as.character(1:(if (is.null(shinyInput$pc)) 1 
                    else ncol(shinyInput$pc))), selected="2"),
                checkboxInput("colbybatchPCA", 
                    "Color By Batch (Default: Color By Condition)", FALSE),
                radioButtons('batchPCA', 'Batch Adjustment',
                    c('None'=0, 'Combat'=1,'SVA'=2), 0),
                width=3
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("PCA", ggvisOutput("plot")),
                    tabPanel("Summary", verbatimTextOutput("PCAsummary")),
                    tabPanel("Table",tableOutput("PCAtable")),
                    tabPanel("Explained Variation",
                        tableOutput("PCAExplainedVariation"))
                ), width=9
            )
        )
    ),
    tabPanel("Shape",
        sidebarLayout(
            sidebarPanel(
                radioButtons('batchShape', 'Batch Adjustment',
                    c('None'=0, 'Combat'=1,'SVA'=2), 0),
                width=3
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Batch Variation", plotOutput("BatchMeanVar"), 
                        br(),
                        h4(paste("Note: Sample-wise p-value is calculated for ",
                        "the variation across samples on the measure",
                        "across genes. Gene-wise p-value is calculated for the",
                        "variation of each gene between batches on the ",
                        "measure across each batch. If the data is quantum ", 
                        "normalized, then the Sample-wise measure across",
                        "genes is same for all samples and Gene-wise p-value",
                        "is a good measure.", 
                        sep=" "))
                    )
                    # tabPanel("Manova", plotOutput("Manova")),
                    # tabPanel("SO Plot", plotOutput("SOplot"))
                ), width=9
            )
        )
    ),
    tabPanel("ComBat",
        sidebarLayout(
            sidebarPanel(
                sliderInput('batches', 'Batch', value=1,
                    min = 1, max = if (is.null(shinyInput$delta.hat)) 1 
                    else nrow(shinyInput$delta.hat), step=1),
                br(),
                radioButtons('optionMeanOnly', 'Batch Adjustment Option',
                    c('Mean and Variance'=0, 'Mean only'=1), 0),
                br(),
                radioButtons('optionParameter', 'Batch Parameter Distribution',
                    c('Parametric'=0, 'Non-Parametric'=1), 0),
                br(),
                actionButton("runCombat", "Run ComBat"),
                shiny::p(paste("Click the button to run ComBat and see the",
                    "results in the other tabs.", sep=" ")),
                width=3
            ),
            mainPanel(
                tabsetPanel(
                    id="CombatMain",
                    tabPanel("ComBat Plots", plotOutput("densityQQPlots", 
                        height = "550px")),
                    tabPanel("Summary", verbatimTextOutput("kstest"), 
                        br(),
                        h4(paste("Note: The non-parametric version of ComBat",
                        "takes much longer time to run and we recommend it",
                        "only when the shape of the non-parametric curve",
                        "widely differs such as a bimodal or highly skewed",
                        "distribution. Otherwise, the difference in batch", 
                        "adjustment is very negligible and parametric version",
                        "is recommended even if p-value of KS test above is",
                        "significant.", sep=" "))
                    ),
                    tabPanel("ComBat Output", value="ComBatOutput", 
                        verbatimTextOutput("combatOutText"))
                ), width=9
            )
        )
    ),
    tabPanel("SVA",
        sidebarLayout(
            sidebarPanel(
                checkboxInput("fsvaOption", 
                    "Frozen SVA (Default: Regression Adjusted)", FALSE),
                actionButton("runSVA", "Run SVA"),
                shiny::p(paste("Click the button to run SVA and see the",
                    "results in the other tabs.", sep=" ")),
                width=3
            ),
            mainPanel(
                tabsetPanel(
                    id="SVAMain",
                    tabPanel("Summary", verbatimTextOutput("svasummary")),
                    tabPanel("SVA Output", value="SVAOutput", 
                        verbatimTextOutput("SVAOutText"))
                ), width=9
            )
        )
    )
))
