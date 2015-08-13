library(shiny)
library(ggvis)
library(d3heatmap)
library(reshape2)

shinyServer(function(input, output, session) {
  PCA <- reactive({
    data.frame(pc[, c(input$xcol,input$ycol)])
  })
  vis_pc <- reactive({
      #interactive PCA plot
      pc$id <- 1:nrow(pc)  
      
      all_values <- function(x) {
        if(is.null(x)) return(NULL)
        row <- pc[pc$id == x$id, ]
        paste0(paste0("PC", input$xcol), ": ", signif(with(row, get(names(row)[input$xcol])),3), "<br>", paste0("PC", input$ycol), ": ", signif(with(row, get(names(row)[input$ycol])),3), "<br>")
      }
    
      pc %>%
        ggvis(~get(names(pc)[input$xcol]), ~get(names(pc)[input$ycol]), fill = ~factor(batch), key := ~id) %>%
        layer_points(size := 75, size.hover := 200) %>%
        add_tooltip(all_values, "hover") %>%
        set_options(hover_duration = 250) %>%
        add_axis("x", title = paste0("PC",input$xcol), properties = axis_props(
          title = list(fontSize = 20),
          labels = list(fontSize = 15)
        )) %>%
        add_axis("y", title = paste0("PC",input$ycol), properties = axis_props(
          title = list(fontSize = 20),
          labels = list(fontSize = 15)
        )) %>%
        add_legend("fill", title = "Batches", properties = legend_props(
          title = list(fontSize = 20),
          labels = list(fontSize = 15)
        ))
  })
  vis_pc %>% bind_shiny("plot")
  output$PCAsummary <- renderPrint({
    summary(pca)
  })
  output$PCAtable <- renderTable({
    PCA()
  })
  
  #interactive boxplot
  vis_bp <- reactive({
    dat <- data.frame(t(data.matrix[1:input$noGenes,]))
    colnames(dat) <- seq(1:input$noGenes)
    dat1 <- melt(dat, id = NULL)
    dat1 %>% ggvis(x = ~variable, y = ~value, fill = ~variable) %>% layer_boxplots() %>%
    add_axis("x", title = paste0("Genes"), properties = axis_props(
      title = list(fontSize = 15),
      labels = list(fontSize = 10)
    )) %>%
      add_axis("y", title = "Expression", title_offset = 75, properties = axis_props(
        title = list(fontSize = 15),
        labels = list(fontSize = 10)
      )) %>%
      add_legend("fill", title = paste0(input$noGenes, " Genes"), properties = legend_props(
        title = list(fontSize = 15),
        labels = list(fontSize = 10)
      ))
  })
  vis_bp %>% bind_shiny("Boxplot")
  data <- reactive({  
    dat <- data.frame(t(data.matrix[1:input$noGenes,]))
    colnames(dat) <- seq(1:input$noGenes)
    dat
  })
  output$BPsummary <- renderPrint({
    summary(data())
  })
  
  output$BPtable <- renderTable({
    data()
  })

  #interactive heatmap
  output$heatmap <- renderD3heatmap({
    d3heatmap(
      data.matrix,
      colors = input$palette,
      labCol = make.unique(as.character(batch)),
      dendrogram = if (input$cluster) "both" else "none"
    ) })
})

#interactive density plots
#   den <- reactive({
#     df <- data.frame(t(gamma.hat))
#     df %>% ggvis(~get(colnames(df)[input$batches])) %>%
#       layer_densities()
#   })
#   den %>% bind_shiny("density")
#})
