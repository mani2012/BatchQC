# Functions to perform surrogate variable analysis within the BatchQC Package 
# Updated 11/30/2015

require(sva)


###### Estimate Surrogate Variables ######

## Determine the number of surrogate variables to estimate in the model using a permutation based procedure
batchQC_num.sv=function(data.matrix, modmatrix){
	data.matrix=as.matrix(data.matrix)
	modmatrix=as.matrix(modmatrix)
	num.sv.be=sva::num.sv(dat=data.matrix,mod=modmatrix, method = "be" , vfilter = NULL, B = 20,seed = 47)
	return(num.sv.be)
}

## Output the data with (known) batch effect removed but biological heterogeneity preserved 
batchQC_psva=function(data.matrix,batch){
	data.matrix=as.matrix(data.matrix)
	batch=as.factor(batch)
	psva.output=psva(dat=data.matrix, batch=batch)
	return(psva.output)
}

## Estimate the surrogate variables using the 2 step approach proposed by Leek and Storey 2007 
batchQC_sva=function(data.matrix, modmatrix){ 
    n.sv=batchQC_num.sv(data.matrix,modmatrix)
    modmatrix0=model.matrix(~1,data=data.frame(t(data.matrix)))
    sva.object=sva::sva(dat=data.matrix, mod=modmatrix, mod0=modmatrix0, n.sv=n.sv)
    return(sva.object)
}


###### Clean Data by Removing Surrogate Variables ######

## Use frozen surrogate variable analysis to remove the surrogate variables inferred from sva 
batchQC_fsva_adjusted=function(data.matrix, modmatrix, sva.object){
	fsva.object=sva::fsva(data.matrix, modmatrix, sva.object, data.matrix,method="fast")
	fsva.adjusted.db=fsva.object$db
	return(fsva.adjusted.db)
}

## Regress the surrogate variables out of the expression data 
batchQC_svregress_adjusted=function(data.matrix,modmatrix,sva.object){
	y=data.matrix
	X=cbind(modmatrix,sva.object$sv)
	Hat=solve(t(X)%*%X)%*%t(X)
 	beta=(Hat%*%t(y))
	P=ncol(modmatrix)
	svregress.adjusted.db=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
	return(svregress.adjusted.db)
}










