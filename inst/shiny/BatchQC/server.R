library(shiny)
library(ggvis)
library(d3heatmap)
library(reshape2)
library(BatchQC)
plotPC <- function(v, d, x, y, ...){
    pcVar <- round((d^2)/sum(d^2)*100,2)
    
    xl <- sprintf(paste("PC ", x, ": %.2f%% variance"), pcVar[x])  
    yl <- sprintf(paste("PC ", y, ": %.2f%% variance"), pcVar[y]) 
    
    plot(v[,x], v[,y], xlab=xl, ylab=yl, ...)
}

shinyServer(function(input, output, session) {
  #needed information from BatchQC
  pc <- shinyInput$pc
  cormat <- shinyInput$cormat
  delta.hat <- shinyInput$delta.hat
  gamma.hat <- shinyInput$gamma.hat
  gamma.bar <- shinyInput$gamma.bar
  lcounts <- shinyInput$lcounts
  a.prior <- shinyInput$a.prior
  b.prior <- shinyInput$b.prior
  batch <- shinyInput$batch
  condition <- shinyInput$condition
  
  setInputs <- function(combatFlag)  {
    if (combatFlag)  {
      shinyInput <<- shinyInputCombat
      pc <<- shinyInput$pc
      cormat <<- shinyInput$cormat
      delta.hat <<- shinyInput$delta.hat
      gamma.hat <<- shinyInput$gamma.hat
      gamma.bar <<- shinyInput$gamma.bar
      lcounts <<- shinyInput$lcounts
      a.prior <<- shinyInput$a.prior
      b.prior <<- shinyInput$b.prior
      batch <<- shinyInput$batch
      condition <<- shinyInput$condition
    } else  {
      shinyInput <<- shinyInputOrig
      pc <<- shinyInput$pc
      cormat <<- shinyInput$cormat
      delta.hat <<- shinyInput$delta.hat
      gamma.hat <<- shinyInput$gamma.hat
      gamma.bar <<- shinyInput$gamma.bar
      lcounts <<- shinyInput$lcounts
      a.prior <<- shinyInput$a.prior
      b.prior <<- shinyInput$b.prior
      batch <<- shinyInput$batch
      condition <<- shinyInput$condition
    }
  }
  #setInputs(FALSE)
  
  #interactive PCA
  PCA <- reactive({
    data.frame(pc[, c(input$xcol,input$ycol)])
  })
  
  #interactive PCA plot
  vis_pc <- reactive({
    if (input$combatPCA)  {
      if (is.null(shinyInputCombat))  {
        session$sendCustomMessage(type = 'testmessage',
                                  message = 'First run ComBat from the ComBat tab')
        updateCheckboxInput(session, "combatPCA", value=FALSE)
      } else  {
        setInputs(TRUE)
      }
    } else  {
      setInputs(FALSE)
    }

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
        add_axis("x", title = paste0("PC",input$xcol), properties = axis_props(
          title = list(fontSize = 15),
          labels = list(fontSize = 10)
        )) %>%
        add_axis("y", title = paste0("PC",input$ycol), properties = axis_props(
          title = list(fontSize = 15),
          labels = list(fontSize = 10)
        )) %>%
        add_legend("fill", title = "Batches", properties = legend_props(
          title = list(fontSize = 15),
          labels = list(fontSize = 10)
        ))
  })
  
  #interactive PCA summary
  vis_pc %>% bind_shiny("plot")
  output$PCAsummary <- renderPrint({
    summary(PCA())
  })
  
  #interactive PCA table
  output$PCAtable <- renderTable({
    PCA()
  })

  #interactive boxplot
  BP <- reactive({
    if (input$combatDE)  {
      if (is.null(shinyInputCombat))  {
        updateCheckboxInput(session, "combatDE", value=FALSE)
      }
    }
    dat <- lcounts
    batch1 <- as.factor(batch)
    batch2 <- split(which(batch == batch1), batch1)
    batch3 <- unlist(lapply(1:length(batch2), function(x) batch2[[x]][1:input$noSamples]))
    dat1 <- dat[,batch3]
    colnames(dat1) <- seq(1:ncol(dat))[batch3]
    dat1
  })
  #interactive Differential Expression boxplot
  DE <- reactive({
    if (input$combatDE)  {
      if (is.null(shinyInputCombat))  {
        updateCheckboxInput(session, "combatDE", value=FALSE)
      }
    }
    dat <- lcounts
    cond1 <- as.factor(condition)
    cond2 <- split(which(condition == cond1), cond1)
    cond3 <- unlist(lapply(1:length(cond2), function(x) cond2[[x]][1:input$ncSamples]))
    dat1 <- dat[,cond3]
    colnames(dat1) <- seq(1:ncol(dat))[cond3]
    dat1
  })
  diffex_bp <- reactive({
    if (input$combatDE)  {
      if (is.null(shinyInputCombat))  {
        session$sendCustomMessage(type = 'testmessage',
                                  message = 'First run ComBat from the ComBat tab')
        updateCheckboxInput(session, "combatDE", value=FALSE)
      } else  {
        setInputs(TRUE)
      }
    } else  {
      setInputs(FALSE)
    }
    if (input$sortbybatch)  {
      batch4 <- split(batch, as.factor(batch))
      batch5 <- unlist(lapply(1:length(batch4), function(x) batch4[[x]][1:input$noSamples]))
      dat1 <- BP()
      dat2 <- melt(as.data.frame(dat1), measure.var=colnames(dat1))
      dat2$batch <- as.factor(unlist(lapply(1:length(batch5), function(x) rep(batch5[x], nrow(dat1)))))
      dat2$condition <- as.factor(unlist(lapply(as.numeric(colnames(dat1)), function(x) rep(condition[x], nrow(dat1)))))
    } else  {
      cond4 <- split(condition, as.factor(condition))
      cond5 <- unlist(lapply(1:length(cond4), function(x) cond4[[x]][1:input$ncSamples]))
      dat1 <- DE()
      dat2 <- melt(as.data.frame(dat1), measure.var=colnames(dat1))
      dat2$condition <- as.factor(unlist(lapply(1:length(cond5), function(x) rep(cond5[x], nrow(dat1)))))
      dat2$batch <- as.factor(unlist(lapply(as.numeric(colnames(dat1)), function(x) rep(batch[x], nrow(dat1)))))
    }
    dat2 %>% group_by(batch) %>%
      ggvis(~variable, ~value, fill = if (input$colbybatch) ~batch else ~condition) %>% layer_boxplots() %>%
      add_tooltip(function(dat2){paste0("Sample: ", dat2$variable, "<br>", 
                                        if (input$colbybatch) "Batch: " else "Condition: ",
                                        if (input$colbybatch) dat2$batch else dat2$condition)}, "hover") %>%
      add_axis("x", title = paste(input$ncSamples, "Sample(s) Per Condition", sep =" "), properties = axis_props(
        title = list(fontSize = 15),
        labels = list(fontSize = 5, angle = 90)
      )) %>%
      add_axis("y", title = "Expression", properties = axis_props(
        title = list(fontSize = 15),
        labels = list(fontSize = 10)
      )) %>%
      add_legend("fill", title = if (input$colbybatch) "Batches" else "Conditions", properties = legend_props(
        title = list(fontSize = 15),
        labels = list(fontSize = 10)
      ))
  })
  diffex_bp %>% bind_shiny("DiffExPlot")
  output$DEsummary <- renderPrint({
    if (input$sortbybatch)  {
      summary(BP())
    } else  {
      summary(DE())
    }
  })
  
  output$DEtable <- renderTable({
    if (input$sortbybatch)  {
      BP()
    } else  {
      DE()
    }
  })
  
  #interactive scatter plot
  output$outliers <- renderPlot({
    BatchQC::batchqc_corscatter(shinyInput$data, shinyInput$batch, mod = shinyInput$mod)
  })

  #interactive heatmap
  output$heatmap <- renderD3heatmap({
    if (input$combatHM)  {
      if (is.null(shinyInputCombat))  {
        session$sendCustomMessage(type = 'testmessage',
                                  message = 'First run ComBat from the ComBat tab')
        updateCheckboxInput(session, "combatHM", value=FALSE)
      } else  {
        setInputs(TRUE)
      }
    } else  {
      setInputs(FALSE)
    }
    d3heatmap(
      lcounts,
      colors = "RdBu",
      labCol = make.unique(as.character(batch)),
      dendrogram = if (input$cluster1) "both" else "none"
    ) })
  
  output$correlation <- renderD3heatmap({
    if (input$combatHM)  {
      if (is.null(shinyInputCombat))  {
        session$sendCustomMessage(type = 'testmessage',
                                  message = 'First run ComBat from the ComBat tab')
        updateCheckboxInput(session, "combatHM", value=FALSE)
      } else  {
        setInputs(TRUE)
      }
    } else  {
      setInputs(FALSE)
    }
    d3heatmap(
      cormat,
      colors = "RdBu",
      labCol = sample,
      labRow = sample,
      dendrogram = if (input$cluster2) "both" else "none"
    ) })
  
  #interactive density plots
  output$densityQQPlots <- renderPlot({
    layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE))
    tmp <- density(gamma.hat[input$batches,])
    plot(tmp,  type='l', main="Density Plot",lwd=2)
    xx <- seq(min(tmp$x), max(tmp$x), length=100)
    lines(xx,dnorm(xx,gamma.bar[input$batches],sqrt(shinyInput$t2[input$batches])), col=2,lwd=2)
    qqnorm(gamma.hat[input$batches,])	
    qqline(gamma.hat[input$batches,], col=2,lwd=2)
    tmp <- density(delta.hat[input$batches,])
    invgam <- 1/rgamma(ncol(delta.hat),a.prior[input$batches],b.prior[input$batches])
    tmp1 <- density(invgam)
    plot(tmp, main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)),lwd=2)
    lines(tmp1, col=2,lwd=2)
    qqplot(delta.hat[input$batches,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')	
    lines(c(0,max(invgam)),c(0,max(invgam)),col=2,lwd=2)	
    title('Q-Q Plot')
  })
  output$kstest <- renderPrint({  
    ks.test(gamma.hat[input$batches,], "pnorm", gamma.bar[input$batches], sqrt(shinyInput$t2[input$batches])) # two-sided, exact
  })
  observe({
    if (input$runCombat > 0)  {
      updateTabsetPanel(session, inputId="CombatMain", selected = "ComBatOutput")
    }
  })

  combatOutText <- eventReactive(input$runCombat, {
    if (is.null(shinyInputCombat))  {
      pdata <- data.frame(batch, condition)
      mod <- model.matrix(~as.factor(shinyInput$condition), data=pdata)
      combat_data <- ComBat(dat=shinyInput$data, batch=shinyInput$batch, mod=mod)
      shinyInput <<- list("data"=combat_data, "batch"=batch, "condition"=condition,
                          "report_dir"=shinyInputOrig$report_dir)
      rmdfile <- system.file("reports/batchqc_report.Rmd", package = "BatchQC")
      report_option_binary="110011111"
      report_option_vector <- unlist(strsplit(as.character(report_option_binary), ""))
      dat <- as.matrix(combat_data)
      outputfile <- rmarkdown::render(rmdfile, output_file="combat_batchqc_report.html", 
                                      output_dir=shinyInput$report_dir)
      shinyInputCombat <<- shinyInput
      "Finished running ComBat. Select the other tabs to view the ComBat results"
    } else  {
      "Already ran ComBat. Select the other tabs to view the ComBat results"
    }
  })
  output$combatOutText <- renderText({
    combatOutText()
  })
  
  myheight <- function(){
    prefix <- 'output_'
    name <- "circos"
    width <- session$clientData[[paste(prefix, name, '_width', sep='')]];
    height <- session$clientData[[paste(prefix, name, '_height', sep='')]];
    if (is.null(width) || is.null(height) || width <= 0 || height <= 0)
      return(NULL)
    if (width<height) {
      return(width)
    } else {
      return(height)
    }
  }
  mywidth <- function(){
    prefix <- 'output_'
    name <- "circos"
    width <- session$clientData[[paste(prefix, name, '_width', sep='')]];
    height <- session$clientData[[paste(prefix, name, '_height', sep='')]];
    if (is.null(width) || is.null(height) || width <= 0 || height <= 0)
      return(NULL)
    if (height<width) {
      return(height)
    } else {
      return(width)
    }
  }
  output$circos <- renderPlot({
    if (input$combatCD)  {
      if (is.null(shinyInputCombat))  {
        session$sendCustomMessage(type = 'testmessage',
                                  message = 'First run ComBat from the ComBat tab')
        updateCheckboxInput(session, "combatCD", value=FALSE)
      } else  {
        setInputs(TRUE)
      }
    } else  {
      setInputs(FALSE)
    }
    my.plot(shinyInput$data, shinyInput$batch, input$AggMethod)
  }, width=mywidth, height=myheight)
})

