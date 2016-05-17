# Simple function to convert binary string to decimal
BinToDec <- function(x) sum(2^(which(rev(unlist(strsplit(as.character(x), 
    "")) == 1)) - 1))

#' Compute log2(counts per mil reads) and library size for each sample
#'
#' @param qcounts quantile normalized counts
#' @param lib.size default is colsums(qcounts)
#' @return list containing log2(quantile counts per mil reads) and library sizes
#' @export
#' @examples
#' nbatch <- 3
#' ncond <- 2
#' npercond <- 10
#' data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, 
#'     npercond=npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, 
#'     seed=1234)
#' data.matrix <- as.matrix(data.matrix)
#' log2CPM(data.matrix)
log2CPM <- function(qcounts, lib.size = NULL) {
    if (is.null(lib.size)) 
        lib.size <- colSums(qcounts)
    minimum <- min(qcounts)
    if (minimum < 0) {
        qcounts <- qcounts - minimum
    }
    avg <- mean(qcounts)
    qcounts <- apply(qcounts, 1:2, FUN = function(x) {
        ifelse(x <= 0, avg, x)
    })
    y <- t(log2(t(qcounts + 0.5)/(lib.size + 1) * 1e+06))
    return(list(y = y, lib.size = lib.size))
}

hcbHeatmap <- function(x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = bluered, colsep, rowsep, sepcolor = "white", 
    sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan", 
    na.color = par("bg"), trace = c("none", "column", "row", 
    "both"), tracecol = "cyan", hline = median(breaks), vline = median(breaks), 
    linecol = tracecol, margins = c(5, 5), ColSideColors, RowSideColors, 
    cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
    labCol = NULL, key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
    ...) {
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none" else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", 
            "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column" else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row" else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    } else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    } else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    } else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    } else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        } else colInd <- rowInd
    } else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    } else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    } else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd] else rownames(x) else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd] else colnames(x) else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    } else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16 else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks) else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.matrix(ColSideColors)) {
                if (is.vector(ColSideColors)) {
                    ColSideColors <- cbind(ColSideColors)
                } else {
                    stop("'ColSideColors' must be a matrix")
                }
            }
            
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc) 
                stop(
                "'ColSideColors' must be a character matrix with ncol(x) rows")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.matrix(RowSideColors)) {
                if (is.vector(RowSideColors)) {
                    RowSideColors <- cbind(RowSideColors)
                } else {
                    stop("'RowSideColors' must be a matrix")
                }
            }
            
            if (!is.character(RowSideColors) || nrow(RowSideColors) != nr) 
                stop(
                "'RowSideColors' must be a character matrix with nrow(x) rows")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), 
                lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        rsc = RowSideColors[rowInd, , drop = FALSE]
        rsc.colors = matrix()
        rsc.names = names(table(rsc))
        rsc.i = 1
        for (rsc.name in rsc.names) {
            rsc.colors[rsc.i] = rsc.name
            rsc[rsc == rsc.name] = rsc.i
            rsc.i = rsc.i + 1
        }
        rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
        image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
        if (length(colnames(RowSideColors)) > 1) {
            axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), 
            colnames(RowSideColors), las = 2, tick = FALSE)
        }
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        csc = ColSideColors[colInd, , drop = FALSE]
        csc.colors = matrix()
        csc.names = names(table(csc))
        csc.i = 1
        for (csc.name in csc.names) {
            csc.colors[csc.i] = csc.name
            csc[csc == csc.name] = csc.i
            csc.i = csc.i + 1
        }
        csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
        image(csc, col = as.vector(csc.colors), axes = FALSE)
        if (length(colnames(ColSideColors)) > 1) {
            axis(2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), 
                colnames(ColSideColors), las = 2, tick = FALSE)
        }
    }
    
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    } else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!is.null(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
            length(csep)), xright = csep + 0.5 + sepwidth[1], 
            ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    } else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    } else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        } else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2) 
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2) 
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        } else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        } else title("Color Key")
    } else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}

batchqc_f.pvalue <- function(dat, mod, mod0) {
    ## F-test (full/reduced model) and returns R2 values
    ## (full/reduced) as well.
    mod00 <- matrix(rep(1, ncol(dat)), ncol = 1)
    n <- dim(dat)[2]
    m <- dim(dat)[1]
    df1 <- dim(mod)[2]
    df0 <- dim(mod0)[2]
    p <- rep(0, m)

    resid <- dat - dat %*% mod %*% solve(t(mod) %*% mod) %*% t(mod)
    rss1 <- rowSums(resid * resid)
    rm(resid)
    
    resid0 <- dat - dat %*% mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0)
    rss0 <- rowSums(resid0 * resid0)
    rm(resid0)

    resid00 <- dat - dat %*% mod00 %*% solve(t(mod00) %*% mod00) %*% t(mod00)
    rss00 <- rowSums(resid00 * resid00)
    rm(resid00)
    
    r2_full <- 1 - rss1/rss00
    r2_reduced <- 1 - rss0/rss00
    
    p <- 1
    if (df1 > df0)  {
        fstats <- ((rss0 - rss1)/(df1 - df0))/(rss1/(n - df1))
        p <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
    }
    return(list(p = p, r2_full = r2_full, r2_reduced = r2_reduced))
}

#' Returns a list of explained variation by batch and condition 
#' combinations
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param condition Condition covariate of interest
#' @param batch Batch covariate 
#' @return List of explained variation by batch and condition
#' @export
#' @examples
#' nbatch <- 3
#' ncond <- 2
#' npercond <- 10
#' data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
#'     npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
#' batch <- rep(1:nbatch, each=ncond*npercond)
#' condition <- rep(rep(1:ncond, each=npercond), nbatch)
#' batchqc_explained_variation(data.matrix, condition, batch)
batchqc_explained_variation <- function(data.matrix, condition, batch) {
    nlb <- nlevels(as.factor(batch))
    nlc <- nlevels(as.factor(condition))
    if ((nlb <= 1)&&(nlc <= 1))  {
        cond_mod <- matrix(rep(1, ncol(data.matrix)), ncol = 1)
        batch_mod <- matrix(rep(1, ncol(data.matrix)), ncol = 1)
    } else if(nlb <= 1)  {
        cond_mod <- model.matrix(~as.factor(condition))
        batch_mod <- matrix(rep(1, ncol(data.matrix)), ncol = 1)
    } else if(nlc <= 1)  {
        cond_mod <- matrix(rep(1, ncol(data.matrix)), ncol = 1)
        batch_mod <- model.matrix(~as.factor(batch))
    } else {
        cond_mod <- model.matrix(~as.factor(condition))
        batch_mod <- model.matrix(~as.factor(batch))
    }
    mod <- cbind(cond_mod, batch_mod[, -1])
    
    cond_test <- batchqc_f.pvalue(data.matrix, mod, batch_mod)
    batch_test <- batchqc_f.pvalue(data.matrix, mod, cond_mod)
    
    cond_ps <- cond_test$p
    batch_ps <- batch_test$p
    
    r2_full <- cond_test$r2_full
    cond_r2 <- batch_test$r2_reduced
    batch_r2 <- cond_test$r2_reduced
    explained_variation <- round(cbind(`Full (Condition+Batch)` = r2_full, 
        Condition = cond_r2, Batch = batch_r2), 5) * 100
    rownames(explained_variation) <- rownames(data.matrix)
    batchqc_ev <- list(explained_variation = explained_variation, 
        cond_test = cond_test, batch_test = batch_test)
    
    return(batchqc_ev)
}

#' Returns explained variation for each principal components
#' 
#' @param pcs Principal components in the given data
#' @param vars Variance of the Principal components in the given data
#' @param condition Condition covariate of interest
#' @param batch Batch covariate 
#' @return Explained variation table for each principal components
#' @export
#' @examples
#' nbatch <- 3
#' ncond <- 2
#' npercond <- 10
#' data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
#'     npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
#' batch <- rep(1:nbatch, each=ncond*npercond)
#' condition <- rep(rep(1:ncond, each=npercond), nbatch)
#' pdata <- data.frame(batch, condition)
#' modmatrix = model.matrix(~as.factor(condition), data=pdata)
#' pca <- batchqc_pca(data.matrix, batch, mod=modmatrix)
#' pcs <- t(data.frame(pca$x))
#' batchqc_pc_explained_variation(pcs, pca$sdev^2, condition, batch)
batchqc_pc_explained_variation <- function(pcs, vars, condition, batch) {
    nlb <- nlevels(as.factor(batch))
    nlc <- nlevels(as.factor(condition))
    if ((nlb <= 1)&&(nlc <= 1))  {
        cond_mod <- matrix(rep(1, ncol(pcs)), ncol = 1)
        batch_mod <- matrix(rep(1, ncol(pcs)), ncol = 1)
    } else if(nlb <= 1)  {
        cond_mod <- model.matrix(~as.factor(condition))
        batch_mod <- matrix(rep(1, ncol(pcs)), ncol = 1)
    } else if(nlc <= 1)  {
        cond_mod <- matrix(rep(1, ncol(pcs)), ncol = 1)
        batch_mod <- model.matrix(~as.factor(batch))
    } else {
        cond_mod <- model.matrix(~as.factor(condition))
        batch_mod <- model.matrix(~as.factor(batch))
    }
    mod <- cbind(cond_mod, batch_mod[, -1])
    cond_test <- batchqc_f.pvalue(pcs, mod, batch_mod)
    batch_test <- batchqc_f.pvalue(pcs, mod, cond_mod)
    cond_ps <- round(cond_test$p, 5)
    batch_ps <- round(batch_test$p, 5)
    r2_full <- round(cond_test$r2_full * 100, 1)
    cond_r2 <- round(batch_test$r2_reduced * 100, 1)
    batch_r2 <- round(cond_test$r2_reduced * 100, 1)
    overlap_r2 <- round((cond_test$r2_reduced + batch_test$r2_reduced - 
        cond_test$r2_full) * 100, 1)
    pvars <- vars*100.0/sum(vars)
    explained_variation <- cbind(`Proportion of Variance (%)` = pvars, 
        `Cumulative Proportion of Variance (%)` = cumsum(pvars), 
        `Percent Variation Explained by Either Condition or Batch` = r2_full, 
        `Percent Variation Explained by Condition` = cond_r2, 
        `Condition Significance (p-value)` = cond_ps, 
        `Percent Variation Explained by Batch` = batch_r2, 
        `Batch Significance (p-value)` = batch_ps)
    rownames(explained_variation) <- rownames(pcs)
    return(explained_variation)
}

#' Returns adjusted data after remove the variation across conditions
#' 
#' @param data.matrix Given data or simulated data from rnaseq_sim()
#' @param batch Batch covariate 
#' @param condition Condition covariate of interest
#' @return Adjusted data after remove the variation across conditions
#' @export
#' @examples
#' nbatch <- 3
#' ncond <- 2
#' npercond <- 10
#' data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=
#'     npercond, ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)
#' batch <- rep(1:nbatch, each=ncond*npercond)
#' condition <- rep(rep(1:ncond, each=npercond), nbatch)
#' batchQC_condition_adjusted(data.matrix, batch, condition)
batchQC_condition_adjusted = function(data.matrix, batch, condition) {
    y <- data.matrix
    nlc <- nlevels(as.factor(condition))
    if (nlc <= 1)  {
        return(Y)
    }
    P <- nlevels(as.factor(batch))
    if (P <= 1)  {
        pdata <- data.frame(condition)
        X <- model.matrix(~as.factor(condition), data = pdata)
        P <- 1
    } else  {
        pdata <- data.frame(batch, condition)
        X <- model.matrix(~as.factor(batch) + as.factor(condition), data=pdata)
    }
    Hat <- solve(t(X) %*% X) %*% t(X)
    beta <- (Hat %*% t(y))
    condition_adjusted <- y - t(as.matrix(X[, -c(1:P)]) %*% beta[-c(1:P), ])
    return(condition_adjusted)
} 
