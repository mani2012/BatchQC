library(shiny)
library(ggvis)
library(d3heatmap)
library(reshape2)
plotPC <- function(v, d, x, y, ...){
    pcVar <- round((d^2)/sum(d^2)*100,2)
    
    xl <- sprintf(paste("PC ", x, ": %.2f%% variance"), pcVar[x])  
    yl <- sprintf(paste("PC ", y, ": %.2f%% variance"), pcVar[y]) 
    
    plot(v[,x], v[,y], xlab=xl, ylab=yl, ...)
}

shinyServer(function(input, output, session) {
  pc <- shinyInput$pc
  cormat <- shinyInput$cormat
  delta.hat <- shinyInput$delta.hat
  gamma.hat <- shinyInput$gamma.hat
  gamma.bar <- shinyInput$gamma.bar
  lcounts <- shinyInput$lcounts
  t2 <- shinyInput$t2
  a.prior <- shinyInput$a.prior
  b.prior <- shinyInput$b.prior
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
        set_options(hover_duration = 150) %>%
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
  
  output$svd <- renderPlot({
    plotPC(res$v,res$d, input$xcol, input$ycol,
           col=cc, # color by condition
           pch=19, main="PCA plot",
           xlim=c(min(res$v[,input$xcol])-.08,max(res$v[,input$xcol])+.08),
           ylim=c(min(res$v[,input$ycol])-.08,max(res$v[,input$ycol])+.08))
    text(res$v[,input$xcol], res$v[,input$ycol], batch, pos=1, cex=0.6)
  }) 
  
  #interactive boxplot
  vis_bp <- reactive({
    dat <- data.frame(data.matrix[,1:input$noSamples])
    colnames(dat) <- seq(1:input$noSamples)
    dat1 <- melt(dat, id = NULL)
    dat1$batch <- as.character(batch)[1:input$noSamples]
    dat1 %>% ggvis(x = ~variable, y = ~value, fill = ~variable) %>% layer_boxplots() %>%
    add_axis("x", title = paste0("Samples"), properties = axis_props(
      title = list(fontSize = 15),
      labels = list(fontSize = 10)
    )) %>%
      add_axis("y", title = "Expression", title_offset = 75, properties = axis_props(
        title = list(fontSize = 15),
        labels = list(fontSize = 10)
      )) %>%
      add_legend("fill", title = paste0(input$noSamples, " Samples"), properties = legend_props(
        title = list(fontSize = 15),
        labels = list(fontSize = 10)
      ))
  })
  vis_bp %>% bind_shiny("Boxplot")
  data <- reactive({  
    dat <- data.frame(data.matrix[,1:input$noSamples])
    colnames(dat) <- seq(1:input$noSamples)
    dat
  })
  output$BPsummary <- renderPrint({
    summary(data())
  })
  
  output$BPtable <- renderTable({
    data()
  })
  
  output$outliers <- renderPlot({
    BatchQC::batchqc_corscatter(data.matrix, batch, mod = shinyInput$mod)
  })

  #interactive heatmap
  output$heatmap <- renderD3heatmap({
    d3heatmap(
      lcounts,
      colors = "RdBu",
      labCol = make.unique(as.character(batch)),
      dendrogram = if (input$cluster1) "both" else "none"
    ) })
  
  output$correlation <- renderD3heatmap({
    d3heatmap(
      cormat,
      colors = "RdBu",
      labCol = sample,
      labRow = sample,
      dendrogram = if (input$cluster2) "both" else "none"
    ) })
  
  output$densityQQPlots <- renderPlot({
    layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE))
    tmp <- density(gamma.hat[input$batches,])
    plot(tmp,  type='l', main="Density Plot")
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    lines(xx,dnorm(xx,gamma.bar[input$batches],sqrt(t2[input$batches])), col=2)
    qqnorm(gamma.hat[input$batches,])	
    qqline(gamma.hat[input$batches,], col=2)
    tmp <- density(delta.hat[input$batches,])
    invgam <- 1/rgamma(ncol(delta.hat),a.prior[input$batches],b.prior[input$batches])
    tmp1 <- density(invgam)
    plot(tmp, main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
    lines(tmp1, col=2)
    qqplot(delta.hat[input$batches,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')	
    lines(c(0,max(invgam)),c(0,max(invgam)),col=2)	
    title('Q-Q Plot')
  })
  output$kstest <- renderPrint({  
    ks.test(gamma.hat[input$batches,], "pnorm", gamma.bar[input$batches], sqrt(t2[input$batches])) # two-sided, exact
  })
})

