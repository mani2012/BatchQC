#  LINEAR MODELS

#' Fit linear model for each gene given a series of arrays. 
#' This is the standard lmFit function from limma package with the modification
#' to accept an additional correlation matrix parameter option
#'
#' @param object A matrix-like data object containing log-ratios or 
#' log-expression values for a series of arrays, with rows corresponding to 
#' genes and columns to samples. Any type of data object that can be processed 
#' by getEAWP is acceptable.
#' @param design the design matrix of the microarray experiment, with rows 
#' corresponding to arrays and columns to coefficients to be estimated. 
#' Defaults to the unit vector meaning that the arrays are treated as replicates
#' @param ndups	positive integer giving the number of times each distinct probe 
#' is printed on each array.
#' @param spacing positive integer giving the spacing between duplicate 
#' occurrences of the same probe, spacing=1 for consecutive rows.
#' @param block	vector or factor specifying a blocking variable on the arrays. 
#' Has length equal to the number of arrays. Must be NULL if ndups>2.
#' @param correlation the inter-duplicate or inter-technical replicate 
#' correlation
#' @param cormatrix the complete correlation matrix of the samples 
#' @param weights non-negative observation weights. Can be a numeric matrix of 
#' individual weights, of same size as the object expression matrix, or a 
#' numeric vector of array weights with length equal to ncol of the expression 
#' matrix, or a numeric vector of gene weights with length equal to nrow of the 
#' expression matrix.
#' @param method fitting method; "ls" for least squares or "robust" for robust 
#' regression
#' @param ... other optional arguments to be passed to lm.series, gls.series or 
#' mrlm
#' @return list containing log2(quantile counts per mil reads) and library sizes
#' @export
lmFitC <- function(object,design=NULL,ndups=1,spacing=1,block=NULL,correlation,
    cormatrix=NULL, weights=NULL,method="ls",...)
#	Fit genewise linear models
#	Gordon Smyth
#	30 June 2003.  Last modified 6 Oct 2015.
{
#	Extract components from y
	y <- getEAWP(object)

#	Check design matrix
	if(is.null(design)) design <- y$design
	if(is.null(design))
		design <- matrix(1,ncol(y$exprs),1)
	else {
		design <- as.matrix(design)
		if(mode(design) != "numeric") stop("design must be a numeric matrix")
		if(nrow(design) != ncol(y$exprs)) stop(
		"row dimension of design doesn't match column dimension of data object")
	}
	ne <- nonEstimable(design)
	if(!is.null(ne)) 
	    cat("Coefficients not estimable:",paste(ne,collapse=" "),"\n")

#	Weights and spacing arguments can be specified in call or stored in y
#	Precedence for these arguments is
#	1. Specified in function call
#	2. Stored in object
#	3. Default values
	if(missing(ndups) && !is.null(y$printer$ndups)) ndups <- y$printer$ndups
	if(missing(spacing) && !is.null(y$printer$spacing)) 
	    spacing <- y$printer$spacing
	if(missing(weights) && !is.null(y$weights)) weights <- y$weights

#	Check method
	method <- match.arg(method,c("ls","robust"))

# If duplicates are present, reduce probe-annotation and Amean to correct length
	if(ndups>1) {
		if(!is.null(y$probes)) 
		    y$probes <- uniquegenelist(y$probes,ndups=ndups,spacing=spacing)
		if(!is.null(y$Amean)) 
		    y$Amean <- rowMeans(unwrapdups(as.matrix(y$Amean),ndups=ndups,
		        spacing=spacing),na.rm=TRUE)
	}

#	Dispatch fitting algorithms
	if(method=="robust")
		fit <- mrlm(y$exprs,design=design,ndups=ndups,spacing=spacing,
	        weights=weights,...)
	else if (!is.null(cormatrix))  {
	    fit <- gls.series.C(y$exprs,design=design,ndups=ndups,spacing=spacing,
            block=block,cormatrix=cormatrix,weights=weights,...)
	} else  {
		if(ndups < 2 && is.null(block))
			fit <- lm.series(y$exprs,design=design,ndups=ndups,spacing=spacing,
			    weights=weights)
		else {
			if(missing(correlation)) 
			    stop("the correlation must be set, see duplicateCorrelation")
			fit <- gls.series(y$exprs,design=design,ndups=ndups,spacing=spacing,
			    block=block,correlation=correlation,weights=weights,...)
		}
	}
	
#	Possible warning on missing coefs
	if(NCOL(fit$coef)>1) {
		n <- rowSums(is.na(fit$coef))
		n <- sum(n>0 & n<NCOL(fit$coef))
		if(n>0) 
		    warning("Partial NA coefficients for ",n," probe(s)",call.=FALSE) 
	}

#	Output
	fit$genes <- y$probes
	fit$Amean <- y$Amean
	fit$method <- method
	fit$design <- design
	new("MArrayLM",fit)
}

gls.series.C <- function(M,design=NULL,ndups=2,spacing=1,block=NULL,
    correlation=NULL,cormatrix=NULL,weights=NULL,...)
#	Fit linear model for each gene to a series of microarrays.
#	Fit is by generalized least squares allowing for correlation between 
#	duplicate spots.
#	Gordon Smyth
#	11 May 2002.  Last revised 29 Dec 2015.
{
#	Check M
	M <- as.matrix(M)
	ngenes <- nrow(M)
	narrays <- ncol(M)

#	Check design
	if(is.null(design)) design <- matrix(1,narrays,1)
	design <- as.matrix(design)
	if(nrow(design) != narrays) 
	    stop("Number of rows of design matrix does not match number of arrays")
	nbeta <- ncol(design)
	coef.names <- colnames(design)

#	Check correlation
	#if(is.null(correlation)) 
	#    correlation <- duplicateCorrelation(M,design=design,ndups=ndups,
	#    spacing=spacing,block=block,weights=weights,...)$consensus.correlation
	#if(abs(correlation) >= 1) 
	#    stop("correlation is 1 or -1, so the model is degenerate")

#	Check weights
	if(!is.null(weights)) {
		weights[is.na(weights)] <- 0
		weights <- asMatrixWeights(weights,dim(M))
		M[weights < 1e-15 ] <- NA
		weights[weights < 1e-15] <- NA
	}

	if(is.null(cormatrix)) {
#	Unwrap duplicates (if necessary) and make correlation matrix
    	if(is.null(block)) {
#		Correlated within-array probes
    		if(ndups<2) {
    			warning("No duplicates (ndups<2)")
    			ndups <- 1
    			correlation <- 0
    		}
    		cormatrix <- diag(rep(correlation,len=narrays),nrow=narrays,
    		    ncol=narrays) %x% array(1,c(ndups,ndups))
    		if(is.null(spacing)) spacing <- 1
    		M <- unwrapdups(M,ndups=ndups,spacing=spacing)
    		if(!is.null(weights)) 
    		    weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
    		design <- design %x% rep(1,ndups)
    		colnames(design) <- coef.names
    		ngenes <- nrow(M)
    		narrays <- ncol(M)
    	} else {
#		Correlated samples
    		ndups <- spacing <- 1
    		block <- as.vector(block)
    		if(length(block)!=narrays) 
    		    stop("Length of block does not match number of arrays")
    		ub <- unique(block)
    		nblocks <- length(ub)
    		Z <- matrix(block,narrays,nblocks)==
    		    matrix(ub,narrays,nblocks,byrow=TRUE)
    		cormatrix <- Z%*%(correlation*t(Z))
    	}
    	diag(cormatrix) <- 1
	}

	stdev.unscaled <- 
	    matrix(NA,ngenes,nbeta,dimnames=list(rownames(M),coef.names))

#	If weights and missing values are absent, use a fast computation
	NoProbeWts <- all(is.finite(M)) && (is.null(weights) || 
	    !is.null(attr(weights,"arrayweights")))
	if(NoProbeWts) {
		V <- cormatrix
		if(!is.null(weights)) {
			wrs <- 1/sqrt(weights[1,])
			V <- wrs * t(wrs * t(V))
		}
		cholV <- chol(V)
		y <- backsolve(cholV,t(M),transpose=TRUE)
		dimnames(y) <- rev(dimnames(M))
		X <- backsolve(cholV,design,transpose=TRUE)
		dimnames(X) <- dimnames(design)
		fit <- lm.fit(X,y)
		if(fit$df.residual>0) {
			if(is.matrix(fit$effects))
				fit$sigma <- sqrt(colMeans(fit$effects[-(1:fit$rank),,
				    drop=FALSE]^2))
			else
				fit$sigma <- sqrt(mean(fit$effects[-(1:fit$rank)]^2))
		} else
			fit$sigma <- rep(NA,ngenes)
		fit$fitted.values <- fit$residuals <- fit$effects <- NULL
		fit$coefficients <- t(fit$coefficients)
		fit$cov.coefficients <- chol2inv(fit$qr$qr,size=fit$qr$rank)
		est <- fit$qr$pivot[1:fit$qr$rank]
		dimnames(fit$cov.coefficients) <- list(coef.names[est],coef.names[est])
		stdev.unscaled[,est] <- matrix(sqrt(diag(fit$cov.coefficients)),ngenes,
		    fit$qr$rank,byrow = TRUE)
		fit$stdev.unscaled <- stdev.unscaled
		fit$df.residual <- rep.int(fit$df.residual,ngenes)
		dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- 
		    dimnames(fit$coefficients)
		fit$pivot <- fit$qr$pivot
		fit$ndups <- ndups
		fit$spacing <- spacing
		fit$block <- block
		fit$correlation <- correlation
		return(fit)
	}

#	Weights or missing values are present, to have to iterate over probes
	beta <- stdev.unscaled
	sigma <- rep(NA,ngenes)
	df.residual <- rep(0,ngenes)
	for (i in 1:ngenes) {
		y <- drop(M[i,])
		o <- is.finite(y)
		y <- y[o]
		n <- length(y)
		if(n > 0) {
			X <- design[o,,drop=FALSE]
			V <- cormatrix[o,o]
			if(!is.null(weights)) {
				wrs <- 1/sqrt(drop(weights[i,o]))
				V <- wrs * t(wrs * t(V))
			}
			cholV <- chol(V)
			y <- backsolve(cholV,y,transpose=TRUE)
			if(all(X==0)) {
				df.residual[i] <- n
				sigma[i] <- sqrt( array(1/n,c(1,n)) %*% y^2 )
			} else {
				X <- backsolve(cholV,X,transpose=TRUE)
				out <- lm.fit(X,y)
				est <- !is.na(out$coefficients)
				beta[i,] <- out$coefficients
				stdev.unscaled[i,est] <- sqrt(diag(chol2inv(out$qr$qr,
				    size=out$rank)))
				df.residual[i] <- out$df.residual
				if(df.residual[i] > 0)
					sigma[i] <- sqrt( array(1/out$df.residual,c(1,n)) %*% 
					    out$residuals^2 )
			}
		}
	}
	cholV <- chol(cormatrix)
	QR <- qr(backsolve(cholV,design,transpose=TRUE))
	cov.coef <- chol2inv(QR$qr,size=QR$rank)
	est <- QR$pivot[1:QR$rank]
	dimnames(cov.coef) <- list(coef.names[est],coef.names[est])
	list(coefficients=beta,stdev.unscaled=stdev.unscaled,sigma=sigma,
	     df.residual=df.residual,ndups=ndups,spacing=spacing,block=block,
	     correlation=correlation,cov.coefficients=cov.coef,pivot=QR$pivot,
	     rank=QR$rank)
}

