library(shiny)
library(pander)
library(ggvis)
library(d3heatmap)
library(reshape2)
library(limma)
library(sva)
# library(HTShape)
library(BatchQC)

plotPC <- function(v, d, x, y, ...) {
    pcVar <- round((d^2)/sum(d^2) * 100, 2)
    
    xl <- sprintf(paste("PC ", x, ": %.2f%% variance"), pcVar[x])
    yl <- sprintf(paste("PC ", y, ": %.2f%% variance"), pcVar[y])
    
    plot(v[, x], v[, y], xlab = xl, ylab = yl, ...)
}

shinyServer(function(input, output, session) {
    # needed information from BatchQC
    shinyInput <- getShinyInput()
    batch <- shinyInput$batch
    condition <- shinyInput$condition
    
    setInputs <- function(batchAdjust) {
        if (batchAdjust == 1) {
            setShinyInput(getShinyInputCombat())
        } else if (batchAdjust == 2) {
            setShinyInput(getShinyInputSVA())
        } else {
            setShinyInput(getShinyInputOrig())
        }
    }
    
    # Summary Table
    catbatch <- function(x) {
        paste("Batch", x)
    }
    catcondition <- function(x) {
        paste("Condition", x)
    }
    output$CondBatchTable <- renderTable({
        counts = table(condition, batch)
        countsmatrix <- as.matrix(counts)
        colnames(countsmatrix) <- sapply(colnames(countsmatrix), catbatch, 
            simplify = TRUE)
        rownames(countsmatrix) <- sapply(rownames(countsmatrix), catcondition, 
            simplify = TRUE)
        countsmatrix
    })
    output$ConfoundTable <- renderTable({
        counts = table(condition, batch)
        rowsums = apply(counts, 1, sum)
        colsums = apply(counts, 2, sum)
        tablesum = sum(rowsums)
        expected = matrix(0, nrow(counts), ncol(counts))
        for (i in 1:nrow(counts)) {
            for (j in 1:ncol(counts)) {
                expected[i, j] = rowsums[i] * colsums[j]/tablesum
            }
        }
        chi = sum((counts - expected)^2/expected)
        mmin = min(nrow(counts), ncol(counts))
        confound1 = sqrt(chi * mmin/((chi + tablesum) * (mmin - 1)))  
        ## Standardized Pearson Correlation Coefficient
        confound2 = sqrt(chi/(tablesum * (mmin - 1)))  ## Cramer's V
        confound <- matrix(c(confound1, confound2), nrow = 1)
        colnames(confound) <- c("Standardized Pearson Correlation Coefficient", 
            "Cramer's V")
        rownames(confound) <- c(
        "Confounding Coefficients (0=no confounding, 1=complete confounding)")
        confound
    })
    
    # Variation Analysis
    output$VariationTable <- renderTable({
        if (input$batchVA == 1) {
            if (is.null(getShinyInputCombat())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run ComBat from the ComBat tab")
                updateRadioButtons(session, "batchVA", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(1)
            }
        } else if (input$batchVA == 2) {
            if (is.null(getShinyInputSVA())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run SVA from the SVA tab")
                updateRadioButtons(session, "batchVA", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        shinyInput <- getShinyInput()
        batchqc_ev <- batchqc_explained_variation(shinyInput$data, condition, 
            batch)
        batchqc_ev$explained_variation[1:input$noGenesVA,]
    })
    # Variation plots
    output$VariationPlot <- renderPlot({
        if (input$batchVA == 1) {
            if (is.null(getShinyInputCombat())) {
            } else {
                setInputs(1)
            }
        } else if (input$batchVA == 2) {
            if (is.null(getShinyInputSVA())) {
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        shinyInput <- getShinyInput()
        batchqc_ev <- batchqc_explained_variation(shinyInput$data, condition, 
            batch)
        apply(batchqc_ev$explained_variation, 2, summary)
        boxplot(batchqc_ev$explained_variation, ylab = 
            "Percent Explained Variation", main = 
            "Percent of Variation Explained by Source")
    })
    
    # P-Value Analysis
    output$PvalueTable <- renderTable({
        if (input$batchPA == 1) {
            if (is.null(getShinyInputCombat())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run ComBat from the ComBat tab")
                updateRadioButtons(session, "batchPA", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(1)
            }
        } else if (input$batchPA == 2) {
            if (is.null(getShinyInputSVA())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run SVA from the SVA tab")
                updateRadioButtons(session, "batchPA", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        shinyInput <- getShinyInput()
        batchqc_ev <- batchqc_explained_variation(shinyInput$data, condition, 
            batch)
        cond_ps <- batchqc_ev$cond_test$p
        batch_ps <- batchqc_ev$batch_test$p
        pvalue_table <- rbind(
            `Batch P-values` = c(summary(batch_ps), 
            `Ps<0.05` = mean(batch_ps <= 0.05)),
            `Condition P-values` = c(summary(cond_ps), 
            `Ps<0.05` = mean(cond_ps <= 0.05)))
        pvalue_table
    })
    # P-Value plots
    output$BatchPvaluePlot <- renderPlot({
        if (input$batchPA == 1) {
            if (is.null(getShinyInputCombat())) {
            } else {
                setInputs(1)
            }
        } else if (input$batchPA == 2) {
            if (is.null(getShinyInputSVA())) {
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        shinyInput <- getShinyInput()
        batchqc_ev <- batchqc_explained_variation(shinyInput$data, condition, 
            batch)
        cond_ps <- batchqc_ev$cond_test$p
        batch_ps <- batchqc_ev$batch_test$p
        nf <- layout(mat = matrix(c(1, 2), 2, 1, byrow = TRUE), 
            height = c(1, 3))
        # par(mar=c(3.1, 3.1, 1.1, 2.1))
        par(mar = c(3, 3, 1, 2))
        boxplot(batch_ps, horizontal = TRUE, outline = TRUE, ylim = c(0, 1), 
            frame = F, col = "green1")
        hist(batch_ps, xlim = c(0, 1), col = "pink", main = "")
        title("Distribution of Batch Effect p-values Across Genes")
    })
    output$ConditionPvaluePlot <- renderPlot({
        if (input$batchPA == 1) {
            if (is.null(getShinyInputCombat())) {
            } else {
                setInputs(1)
            }
        } else if (input$batchPA == 2) {
            if (is.null(getShinyInputSVA())) {
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        shinyInput <- getShinyInput()
        batchqc_ev <- batchqc_explained_variation(shinyInput$data, condition, 
            batch)
        cond_ps <- batchqc_ev$cond_test$p
        nf <- layout(mat = matrix(c(1, 2), 2, 1, byrow = TRUE), 
            height = c(1, 3))
        # par(mar=c(3.1, 3.1, 1.1, 2.1))
        par(mar = c(3, 3, 1, 2))
        boxplot(cond_ps, horizontal = TRUE, outline = TRUE, ylim = c(0, 1), 
            frame = F, col = "green1")
        hist(cond_ps, xlim = c(0, 1), col = "pink", main = "")
        title("Distribution of Condition Effect p-values Across Genes")
    })
    
    
    # interactive PCA
    PCA <- reactive({
        if (input$batchPCA == 1) {
            if (is.null(getShinyInputCombat())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run ComBat from the ComBat tab")
                updateRadioButtons(session, "batchPCA", choices = 
                    list(None = 0, Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(1)
            }
        } else if (input$batchPCA == 2) {
            if (is.null(getShinyInputSVA())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run SVA from the SVA tab")
                updateRadioButtons(session, "batchPCA", choices = 
                    list(None = 0, Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        shinyInput <- getShinyInput()
        pc <- shinyInput$pc
        data.frame(pc[, c(input$xcol, input$ycol)])
    })
    
    # interactive PCA plot
    vis_pc <- reactive({
        if (input$batchPCA == 1) {
            if (is.null(getShinyInputCombat())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run ComBat from the ComBat tab")
                updateRadioButtons(session, "batchPCA", choices = 
                    list(None = 0, Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(1)
            }
        } else if (input$batchPCA == 2) {
            if (is.null(getShinyInputSVA())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run SVA from the SVA tab")
                updateRadioButtons(session, "batchPCA", choices = 
                    list(None = 0, Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        shinyInput <- getShinyInput()
        pc <- shinyInput$pc
        
        pc$id <- 1:nrow(pc)
        
        all_values <- function(x) {
            if (is.null(x)) 
                return(NULL)
            row <- pc[pc$id == x$id, ]
            paste0(paste0("PC", input$xcol), ": ", signif(with(row, 
                get(names(row)[input$xcol])), 3), "<br>", 
                paste0("PC", input$ycol), ": ", signif(with(row, 
                get(names(row)[input$ycol])), 3), "<br>")
        }
        
        pc %>% 
        ggvis(~get(names(pc)[input$xcol]), ~get(names(pc)[input$ycol]), 
            fill = if (input$colbybatchPCA) ~factor(batch) 
            else ~factor(condition), `:=`(key, ~id)) %>% 
        layer_points(`:=`(size, 75), `:=`(size.hover, 200)) %>% 
        add_tooltip(all_values, "hover") %>% 
        add_axis("x", title = paste0("PC", input$xcol), properties = axis_props(
            title = list(fontSize = 15), labels = list(fontSize = 10))) %>% 
        add_axis("y", title = paste0("PC", input$ycol), properties = axis_props(
            title = list(fontSize = 15), labels = list(fontSize = 10))) %>% 
        add_legend("fill", title = if (input$colbybatchPCA) 
            "Batches" else "Conditions", properties = legend_props(title = 
            list(fontSize = 15), labels = list(fontSize = 10)))
    })
    
    # interactive PCA summary
    vis_pc %>% bind_shiny("plot")
    output$PCAsummary <- renderPrint({
        summary(PCA())
    })
    
    # interactive PCA table
    output$PCAtable <- renderTable({
        PCA()
    })
    
    output$PCAExplainedVariation <- renderTable({
        PCA()
        shinyInput <- getShinyInput()
        pc <- shinyInput$pc
        pcs <- t(pc)
        explained_variation <- batchqc_pc_explained_variation(pcs, 
            shinyInput$vars, condition, batch)
        explained_variation
    })
    
    # interactive boxplot
    BP <- reactive({
        if (input$batchDE == 1) {
            if (is.null(getShinyInputCombat())) {
                updateRadioButtons(session, "batchDE", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            }
        } else if (input$batchDE == 2) {
            if (is.null(getShinyInputSVA())) {
                updateRadioButtons(session, "batchDE", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            }
        }
        shinyInput <- getShinyInput()
        lcounts <- shinyInput$lcounts
        dat <- lcounts
        batch1 <- as.factor(batch)
        batch2 <- split(which(batch == batch1), batch1)
        batch3 <- unlist(lapply(1:length(batch2), 
            function(x) batch2[[x]][1:input$noSamples]))
        dat1 <- dat[, batch3]
        colnames(dat1) <- seq(1:ncol(dat))[batch3]
        dat1
    })
    # interactive Differential Expression boxplot
    DE <- reactive({
        if (input$batchDE == 1) {
            if (is.null(getShinyInputCombat())) {
                updateRadioButtons(session, "batchDE", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            }
        } else if (input$batchDE == 2) {
            if (is.null(getShinyInputSVA())) {
                updateRadioButtons(session, "batchDE", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            }
        }
        shinyInput <- getShinyInput()
        lcounts <- shinyInput$lcounts
        dat <- lcounts
        cond1 <- as.factor(condition)
        cond2 <- split(which(condition == cond1), cond1)
        cond3 <- unlist(lapply(1:length(cond2), 
            function(x) cond2[[x]][1:input$ncSamples]))
        dat1 <- dat[, cond3]
        colnames(dat1) <- seq(1:ncol(dat))[cond3]
        dat1
    })
    diffex_bp <- reactive({
        if (input$batchDE == 1) {
            if (is.null(getShinyInputCombat())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run ComBat from the ComBat tab")
                updateRadioButtons(session, "batchDE", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(1)
            }
        } else if (input$batchDE == 2) {
            if (is.null(getShinyInputSVA())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run SVA from the SVA tab")
                updateRadioButtons(session, "batchDE", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        if (input$sortbybatch) {
            batch4 <- split(batch, as.factor(batch))
            batch5 <- unlist(lapply(1:length(batch4), 
                function(x) batch4[[x]][1:input$noSamples]))
            dat1 <- BP()
            dat2 <- melt(as.data.frame(dat1), measure.var = colnames(dat1))
            dat2$batch <- as.factor(unlist(lapply(1:length(batch5), 
                function(x) rep(batch5[x], nrow(dat1)))))
            dat2$condition <- as.factor(unlist(lapply(as.numeric(colnames(dat1))
                , function(x) rep(condition[x], nrow(dat1)))))
            dat2$samples <- unlist(lapply(seq(ncol(dat1)), 
                function(x) rep(x, nrow(dat1))))
        } else {
            cond4 <- split(condition, as.factor(condition))
            cond5 <- unlist(lapply(1:length(cond4), 
                function(x) cond4[[x]][1:input$ncSamples]))
            dat1 <- DE()
            dat2 <- melt(as.data.frame(dat1), measure.var = colnames(dat1))
            dat2$condition <- as.factor(unlist(lapply(1:length(cond5), 
                function(x) rep(cond5[x], nrow(dat1)))))
            dat2$batch <- as.factor(unlist(lapply(as.numeric(colnames(dat1)), 
                function(x) rep(batch[x], nrow(dat1)))))
            dat2$samples <- unlist(lapply(seq(ncol(dat1)), 
                function(x) rep(x, nrow(dat1))))
        }
        dat2 %>% group_by(batch) %>% ggvis(~samples, ~value, fill = 
            if (input$colbybatch) ~batch else ~condition) %>%
            layer_boxplots() %>% 
            add_tooltip(function(dat2) { paste0("Sample: ", 
                if (input$sortbybatch) colnames(shinyInput$lcounts)[
                as.numeric(colnames(BP()))[dat2$samples]] else colnames(
                shinyInput$lcounts)[as.numeric(colnames(DE()))[dat2$samples]],
                "<br>", if (input$colbybatch) "Batch: " else "Condition: ", 
                if (input$colbybatch) dat2$batch else dat2$condition)
                }, "hover") %>% 
            add_axis("x", title = if (input$sortbybatch) 
                paste(input$noSamples, "Sample(s) Per Batch", sep = " ") 
                else paste(input$ncSamples, "Sample(s) Per Condition", sep=" "), 
                properties = axis_props(title = list(fontSize = 15), 
                labels = list(fontSize = 5, angle = 90))) %>% 
            add_axis("y", title = "Expression", properties = axis_props(title = 
                list(fontSize = 15),labels = list(fontSize = 10))) %>% 
            add_legend("fill", title = if (input$colbybatch) 
                "Batches" else "Conditions", properties = legend_props(title = 
                list(fontSize = 15), labels = list(fontSize = 10)))
    })
    diffex_bp %>% bind_shiny("DiffExPlot")
    output$DEsummary <- renderPrint({
        if (input$batchDE == 1) {
            if (is.null(getShinyInputCombat())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run ComBat from the ComBat tab")
                updateRadioButtons(session, "batchDE", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(1)
            }
        } else if (input$batchDE == 2) {
            if (is.null(getShinyInputSVA())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run SVA from the SVA tab")
                updateRadioButtons(session, "batchDE", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        if (input$sortbybatch) {
            dat1 <- BP()
            colnames(dat1) <- colnames(shinyInput$lcounts)[as.numeric(
                colnames(dat1))]
            summary(dat1)
        } else {
            dat1 <- DE()
            colnames(dat1) <- colnames(shinyInput$lcounts)[as.numeric(
                colnames(dat1))]
            summary(dat1)
        }
    }, width=80)
    
    # output$DEtable <- renderTable({
    #     if (input$batchDE == 1) {
    #         if (is.null(getShinyInputCombat())) {
    #             session$sendCustomMessage(type = "testmessage", message = 
    #                 "First run ComBat from the ComBat tab")
    #             updateRadioButtons(session, "batchDE", choices = list(None = 0, 
    #                 Combat = 1, SVA = 2), selected = 0)
    #         } else {
    #             setInputs(1)
    #         }
    #     } else if (input$batchDE == 2) {
    #         if (is.null(getShinyInputSVA())) {
    #             session$sendCustomMessage(type = "testmessage", message = 
    #                 "First run SVA from the SVA tab")
    #             updateRadioButtons(session, "batchDE", choices = list(None = 0, 
    #                 Combat = 1, SVA = 2), selected = 0)
    #         } else {
    #             setInputs(2)
    #         }
    #     } else {
    #         setInputs(0)
    #     }
    #     if (input$sortbybatch) {
    #         BP()
    #     } else {
    #         DE()
    #     }
    # })
    
    output$LimmaTable <- renderTable({
        if (input$batchDE == 1) {
            if (is.null(getShinyInputCombat())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run ComBat from the ComBat tab")
                updateRadioButtons(session, "batchDE", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(1)
            }
        } else if (input$batchDE == 2) {
            if (is.null(getShinyInputSVA())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run SVA from the SVA tab")
                updateRadioButtons(session, "batchDE", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        shinyInput <- getShinyInput()
        pdata <- data.frame(batch, condition)
        ncond <- nlevels(as.factor(condition))
        nbatch <- nlevels(as.factor(batch))
        limmaTable <- NULL
        if (ncond > 1)  {
            if (nbatch <= 1)  {
                mod_full <- model.matrix(~as.factor(condition), data = pdata)
            } else  {
                mod_full <- model.matrix(~as.factor(condition) + 
                    ~as.factor(batch), data = pdata)
            }
            fit <- lmFit(shinyInput$data, mod)
            fit2 <- eBayes(fit)
            ncond <- nlevels(as.factor(condition))
            limmaTable <- topTable(fit2, coef = 2:ncond, number = input$noGenes)
            for (j in 2:ncond)  {
                colnames(limmaTable)[j-1] <- paste("Condition: ", 
                levels(as.factor(condition))[j], " (logFC)", sep='')
            }
        }
        limmaTable
    })
    
    # interactive scatter plot
    output$outliers <- renderPlot({
        shinyInput <- getShinyInput()
        BatchQC::batchqc_corscatter(shinyInput$data, shinyInput$batch, 
            mod = shinyInput$mod)
    })
    
    # interactive heatmap
    output$heatmap <- renderD3heatmap({
        if (input$batchHM1 == 1) {
            if (is.null(getShinyInputCombat())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run ComBat from the ComBat tab")
                updateRadioButtons(session, "batchHM1", choices = list(None = 0,
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(1)
            }
        } else if (input$batchHM1 == 2) {
            if (is.null(getShinyInputSVA())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run SVA from the SVA tab")
                updateRadioButtons(session, "batchHM1", choices = list(None = 0,
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        shinyInput <- getShinyInput()
        lcounts <- shinyInput$lcounts
        lcountsReduced <- lcounts[1:input$noGenesHM,]
        batch <- shinyInput$batch
        bc <- rainbow(max(batch))
        colorfun <- function(i) {
            return(bc[i])
        }
        cc <- sapply(batch, colorfun, simplify = TRUE)
        d3heatmap(lcountsReduced, colors = "RdBu", labCol = make.unique(
            as.character(batch)), dendrogram = 
            if (input$cluster1) "both" else "none",
            ColSideColors=cc)
    })
    
    output$correlation <- renderD3heatmap({
        if (input$batchHM2 == 1) {
            if (is.null(getShinyInputCombat())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run ComBat from the ComBat tab")
                updateRadioButtons(session, "batchHM2", choices = list(None = 0,
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(1)
            }
        } else if (input$batchHM2 == 2) {
            if (is.null(getShinyInputSVA())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run SVA from the SVA tab")
                updateRadioButtons(session, "batchHM2", choices = list(None = 0,
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        shinyInput <- getShinyInput()
        cormat <- shinyInput$cormat
        nsample <- dim(shinyInput$data)[2]
        sample <- 1:nsample
        fbatch <- as.factor(shinyInput$batch)
        nbatch <- nlevels(fbatch)
        bc <- rainbow(nbatch)
        intbatch <- as.integer(fbatch)
        colorfun <- function(i) {
            return(bc[i])
        }
        cc <- sapply(intbatch, colorfun, simplify = TRUE)
        d3heatmap(cormat, colors = "RdBu", labCol = sample, labRow = sample, 
            dendrogram = if (input$cluster2) "both" else "none",
            ColSideColors=cc)
    })
    
    # Shape plots
    # output$SOplot <- renderPlot({
    #     if (input$batchShape == 1) {
    #         if (is.null(getShinyInputCombat())) {
    #             session$sendCustomMessage(type = "testmessage", message = 
    #                 "First run ComBat from the ComBat tab")
    #             updateRadioButtons(session, "batchShape", choices = list(
    #                 None = 0, Combat = 1, SVA = 2), selected = 0)
    #         } else {
    #             setInputs(1)
    #         }
    #     } else if (input$batchShape == 2) {
    #         if (is.null(getShinyInputSVA())) {
    #             session$sendCustomMessage(type = "testmessage", message = 
    #                 "First run SVA from the SVA tab")
    #             updateRadioButtons(session, "batchShape", choices = list(
    #                 None = 0, Combat = 1, SVA = 2), selected = 0)
    #         } else {
    #             setInputs(2)
    #         }
    #     } else {
    #         setInputs(0)
    #     }
    #     shinyInput <- getShinyInput()
    #     lcounts <- shinyInput$lcounts
    #     lcounts_adj <- batchQC_condition_adjusted(lcounts, batch, condition)
    #     res <- HTShape::fitShape(lcounts_adj, nLmom = 4)
    #     # t3 <- res$lrats[, 't3'] # Grab L-skew estimates. 
    #     # t4 <- res$lrats[, 't4'] # Grab L-kurt estimates.
    #     t3 <- res$lrats["LR3", ]  # Grab L-skew estimates.
    #     t4 <- res$lrats["LR4", ]  # Grab L-kurt estimates.
    #     HTShape::plotSO(t3, t4, verbose = FALSE)
    # })
    # output$Manova <- renderPlot({
    #     if (input$batchShape == 1) {
    #         if (is.null(getShinyInputCombat())) {
    #             session$sendCustomMessage(type = "testmessage", message = 
    #                 "First run ComBat from the ComBat tab")
    #             updateRadioButtons(session, "batchShape", choices = list(
    #                 None = 0, Combat = 1, SVA = 2), selected = 0)
    #         } else {
    #             setInputs(1)
    #         }
    #     } else if (input$batchShape == 2) {
    #         if (is.null(getShinyInputSVA())) {
    #             session$sendCustomMessage(type = "testmessage", message = 
    #                 "First run SVA from the SVA tab")
    #             updateRadioButtons(session, "batchShape", choices = list(
    #                 None = 0, Combat = 1, SVA = 2), selected = 0)
    #         } else {
    #             setInputs(2)
    #         }
    #     } else {
    #         setInputs(0)
    #     }
    #     shinyInput <- getShinyInput()
    #     lcounts <- shinyInput$lcounts
    #     lcounts_adj <- batchQC_condition_adjusted(lcounts, batch, condition)
    #     bf <- as.factor(shinyInput$batch)
    #     HTShape::shapeManova(lcounts_adj, batch, lrats = TRUE, plot = TRUE, 
    #         groupCol = rainbow(nlevels(bf))[bf])
    # })
    output$BatchMeanVar <- renderPlot({
        if (input$batchShape == 1) {
            if (is.null(getShinyInputCombat())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run ComBat from the ComBat tab")
                updateRadioButtons(session, "batchShape", choices = list(
                    None = 0, Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(1)
            }
        } else if (input$batchShape == 2) {
            if (is.null(getShinyInputSVA())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run SVA from the SVA tab")
                updateRadioButtons(session, "batchShape", choices = list(
                    None = 0, Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        shinyInput <- getShinyInput()
        lcounts <- shinyInput$lcounts
        lcounts_adj <- batchQC_condition_adjusted(lcounts, batch, condition)
        bf <- as.factor(shinyInput$batch)
        batchQC_shapeVariation(lcounts_adj, batch, plot = TRUE, groupCol = 
            rainbow(nlevels(bf))[bf])
    })
    
    # interactive density plots
    output$densityQQPlots <- renderPlot({
        shinyInput <- getShinyInputOrig()
        delta.hat <- shinyInput$delta.hat
        if (!is.null(delta.hat))  {
            gamma.hat <- shinyInput$gamma.hat
            gamma.bar <- shinyInput$gamma.bar
            a.prior <- shinyInput$a.prior
            b.prior <- shinyInput$b.prior
            layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
            tmp <- density(gamma.hat[input$batches, ])
            xx <- seq(min(tmp$x), max(tmp$x), length = 100)
            tmp1 <- dnorm(xx, gamma.bar[input$batches], 
                          sqrt(shinyInput$t2[input$batches]))
            plot(tmp, type = "l", main = "Density Plot", 
                 ylim = c(0, max(tmp$y, tmp1)), lwd = 2)
            lines(xx, tmp1, col = 2, lwd = 2)
            qqnorm(gamma.hat[input$batches, ])
            qqline(gamma.hat[input$batches, ], col = 2, lwd = 2)
            tmp <- density(delta.hat[input$batches, ])
            invgam <- 1/rgamma(ncol(delta.hat), a.prior[input$batches], 
                               b.prior[input$batches])
            tmp1 <- density(invgam)
            plot(tmp, main = "Density Plot", ylim = c(0, max(tmp$y, tmp1$y)), 
                 lwd = 2)
            lines(tmp1, col = 2, lwd = 2)
            qqplot(delta.hat[input$batches, ], invgam, xlab="Sample Quantiles", 
                   ylab = "Theoretical Quantiles")
            lines(c(0, max(invgam)), c(0, max(invgam)), col = 2, lwd = 2)
            title("Q-Q Plot")
        }
    })
    output$kstest <- renderPrint({
        shinyInput <- getShinyInputOrig()
        delta.hat <- shinyInput$delta.hat
        summarytext <- ""
        if (!is.null(delta.hat))  {
            gamma.hat <- shinyInput$gamma.hat
            gamma.bar <- shinyInput$gamma.bar
            a.prior <- shinyInput$a.prior
            b.prior <- shinyInput$b.prior
            ksout <- ks.test(gamma.hat[input$batches, ], "pnorm", 
                             gamma.bar[input$batches], sqrt(shinyInput$t2[input$batches]))  
            # two-sided, exact
            summarytext <- 
                "Batch mean distribution across genes: Normal vs Empirical distribution"
            summarytext <- paste(summarytext, "Two-sided Kolmogorov-Smirnov test", 
                                 sep = "\n")
            summarytext <- paste(summarytext, "Selected Batch: ", sep = "\n")
            summarytext <- paste(summarytext, input$batches, sep = "")
            summarytext <- paste(summarytext, "Statistic D = ", sep = "\n")
            summarytext <- paste(summarytext, signif(ksout$statistic, 4), sep = "")
            summarytext <- paste(summarytext, "p-value = ", sep = "\n")
            summarytext <- paste(summarytext, signif(ksout$p.value, 4), sep = "")
            
            invgam <- 1/rgamma(ncol(delta.hat), a.prior[input$batches], 
                               b.prior[input$batches])
            ksvarout <- ks.test(delta.hat[input$batches, ], invgam)  
            # two-sided, exact
            summarytext <- paste(summarytext, 
                "\n\n\nBatch Variance distribution across genes: ",
                "Inverse Gamma vs Empirical distribution", 
                sep = "")
            summarytext <- paste(summarytext, 
                "Two-sided Kolmogorov-Smirnov test", sep = "\n")
            summarytext <- paste(summarytext, "Selected Batch: ", sep = "\n")
            summarytext <- paste(summarytext, input$batches, sep = "")
            summarytext <- paste(summarytext, "Statistic D = ", sep = "\n")
            summarytext <- paste(summarytext, signif(ksvarout$statistic, 4), 
                sep = "")
            summarytext <- paste(summarytext, "p-value = ", sep = "\n")
            summarytext <- paste(summarytext, signif(ksvarout$p.value, 4), 
                sep = "")
        }
        cat(summarytext)
    })
    observe({
        if (input$runCombat > 0) {
            updateTabsetPanel(session, inputId = "CombatMain", 
                selected = "ComBatOutput")
        }
    })
    
    combatOutText <- eventReactive(input$runCombat, {
        if (is.null(getShinyInputCombat())) {
            outText <- paste("Finished running ComBat.",
                "Select the other tabs to view the ComBat results", sep=" ")
        } else {
            outText <- paste("Finished re-running ComBat.",
                "Select the other tabs to view the ComBat results", sep=" ")
        }
        shinyInput <- getShinyInputOrig()
        batch <- shinyInput$batch
        condition <- shinyInput$condition
        pdata <- data.frame(batch, condition)
        mean.only = FALSE
        if (input$optionMeanOnly == 1) {
            mean.only = TRUE
        }
        par.prior = TRUE
        if (input$optionParameter == 1) {
            par.prior = FALSE
        }
        ncond <- nlevels(as.factor(condition))
        if (ncond <= 1)  {
            mod = matrix(rep(1, ncol(shinyInput$data)), ncol = 1)
        } else  {
            mod = model.matrix(~as.factor(condition), data = pdata)
        }
        #combat_data <- ComBat(dat = shinyInput$data, batch = shinyInput$batch, 
        #    mod = mod, par.prior = par.prior, mean.only = mean.only)
        
        lcounts <- shinyInput$lcounts
        lcounts_adj <- batchQC_condition_adjusted(lcounts, batch, condition)
        gnormdata <- gnormalize(lcounts_adj)
        combat_data <- ComBat(dat=lcounts_adj, batch = shinyInput$batch, 
            mod = mod, par.prior = par.prior, mean.only = mean.only)
        
        report_option_binary = "111111111"
        report_option_vector <- unlist(strsplit(as.character(
            report_option_binary), ""))
        shinyInput <- list(data = combat_data, batch = batch, condition = 
            condition, report_dir = shinyInput$report_dir, 
            report_option_vector = report_option_vector)
        setShinyInput(shinyInput)
        rmdfile <- system.file("reports/batchqc_report.Rmd", 
            package = "BatchQC")
        # dat <- as.matrix(combat_data)
        outputfile <- rmarkdown::render(rmdfile, output_file = 
            "combat_batchqc_report.html", output_dir = shinyInput$report_dir)
        shinyInput <- getShinyInput()
        setShinyInputCombat(shinyInput)
        outText
    })
    output$combatOutText <- renderText({
        combatOutText()
    })
    output$svasummary <- renderText({
        shinyInput <- getShinyInputOrig()
        condition <- shinyInput$condition
        nsample <- dim(shinyInput$data)[2]
        sample <- 1:nsample
        pdata <- data.frame(sample, condition)
        ncond <- nlevels(as.factor(condition))
        if (ncond <= 1)  {
            modmatrix = matrix(rep(1, nsample), ncol = 1)
        } else  {
            modmatrix = model.matrix(~as.factor(condition), data = pdata)
        }
        n.sv <- batchQC_num.sv(shinyInput$data, modmatrix)
        paste("Number of Surrogate Variables found in the given data:", n.sv)
    })
    observe({
        if (input$runSVA > 0) {
            updateTabsetPanel(session, inputId = "SVAMain", 
                selected = "SVAOutput")
        }
    })
    SVAOutText <- eventReactive(input$runSVA, {
        shinyInput <- getShinyInputOrig()
        batch <- shinyInput$batch
        condition <- shinyInput$condition
        if (input$fsvaOption) {
            if (is.null(getShinyInputSVAf())) {
                pdata <- data.frame(batch, condition)
                ncond <- nlevels(as.factor(condition))
                if (ncond <= 1)  {
                    mod = matrix(rep(1, ncol(shinyInput$data)), ncol = 1)
                } else  {
                    mod = model.matrix(~as.factor(condition), data = pdata)
                }
                sva.object <- batchQC_sva(shinyInput$data, mod)
                svaf_data <- batchQC_fsva_adjusted(shinyInput$data, mod, 
                    sva.object)
                report_option_binary = "111111111"
                report_option_vector <- unlist(strsplit(as.character(
                    report_option_binary), ""))
                shinyInput <- list(data = svaf_data, batch = batch, condition = 
                    condition, report_dir = shinyInput$report_dir, 
                    report_option_vector = report_option_vector)
                setShinyInput(shinyInput)
                rmdfile <- system.file("reports/batchqc_report.Rmd", 
                    package = "BatchQC")
                # dat <- as.matrix(svaf_data)
                outputfile <- rmarkdown::render(rmdfile, output_file = 
                    "svaf_batchqc_report.html", output_dir = 
                    shinyInput$report_dir)
                shinyInput <- getShinyInput()
                setShinyInput(shinyInput)
                setShinyInputSVAf(shinyInput)
                setShinyInputSVA(shinyInput)
                paste("Finished running SVA with frozen SVA option.",
                    "Select the other tabs to view the SVA results.", sep=" ")
            } else {
                setShinyInputSVA(getShinyInputSVAf())
                paste("Already ran SVA with frozen SVA option.",
                    "Select the other tabs to view the SVA results.", sep=" ")
            }
        } else {
            if (is.null(getShinyInputSVAr())) {
                pdata <- data.frame(batch, condition)
                ncond <- nlevels(as.factor(condition))
                if (ncond <= 1)  {
                    mod = matrix(rep(1, ncol(shinyInput$data)), ncol = 1)
                } else  {
                    mod = model.matrix(~as.factor(condition), data = pdata)
                }
                sva.object <- batchQC_sva(shinyInput$data, mod)
                svar_data <- batchQC_svregress_adjusted(shinyInput$data, mod, 
                    sva.object)
                report_option_binary = "111111111"
                report_option_vector <- unlist(strsplit(as.character(
                    report_option_binary), ""))
                shinyInput <- list(data = svar_data, batch = batch, condition = 
                    condition, report_dir = shinyInput$report_dir, 
                    report_option_vector = report_option_vector)
                setShinyInput(shinyInput)
                rmdfile <- system.file("reports/batchqc_report.Rmd", 
                    package = "BatchQC")
                # dat <- as.matrix(svar_data)
                outputfile <- rmarkdown::render(rmdfile, output_file = 
                    "svar_batchqc_report.html", output_dir = 
                    shinyInput$report_dir)
                shinyInput <- getShinyInput()
                setShinyInput(shinyInput)
                setShinyInputSVAr(shinyInput)
                setShinyInputSVA(shinyInput)
                paste("Finished running SVA with regression SVA option.",
                    "Select the other tabs to view the SVA results.", sep=" ")
            } else {
                setShinyInputSVA(getShinyInputSVAr())
                paste("Already ran SVA with regression SVA option.",
                    "Select the other tabs to view the SVA results.", sep=" ")
            }
        }
    })
    output$SVAOutText <- renderText({
        SVAOutText()
    })
    
    myheight <- function() {
        prefix <- "output_"
        name <- "circos"
        width <- session$clientData[[paste(prefix, name, "_width", sep = "")]]
        height <- session$clientData[[paste(prefix, name, "_height", sep = "")]]
        if (is.null(width) || is.null(height) || width <= 0 || height <= 0) 
            return(NULL)
        if (width < height) {
            return(width)
        } else {
            return(height)
        }
    }
    mywidth <- function() {
        prefix <- "output_"
        name <- "circos"
        width <- session$clientData[[paste(prefix, name, "_width", sep = "")]]
        height <- session$clientData[[paste(prefix, name, "_height", sep = "")]]
        if (is.null(width) || is.null(height) || width <= 0 || height <= 0) 
            return(NULL)
        if (height < width) {
            return(height)
        } else {
            return(width)
        }
    }
    output$circos <- renderPlot({
        if (input$batchCD == 1) {
            if (is.null(getShinyInputCombat())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run ComBat from the ComBat tab")
                updateRadioButtons(session, "batchCD", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(1)
            }
        } else if (input$batchCD == 2) {
            if (is.null(getShinyInputSVA())) {
                session$sendCustomMessage(type = "testmessage", message = 
                    "First run SVA from the SVA tab")
                updateRadioButtons(session, "batchCD", choices = list(None = 0, 
                    Combat = 1, SVA = 2), selected = 0)
            } else {
                setInputs(2)
            }
        } else {
            setInputs(0)
        }
        shinyInput <- getShinyInput()
        batchqc_circosplot(shinyInput$data, if (input$colbybatchCD) 
            shinyInput$batch else shinyInput$condition, input$AggMethod)
    }, width = mywidth, height = myheight)
}) 
