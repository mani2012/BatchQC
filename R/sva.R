# Determine the number of surrogate variables to estimate in the model
    ## By default uses the permutation based procedure (method="be")
# Estimate the surrogate variables using the 2 step approach proposed by Leek and Storey 2007 
	## psva: output data with batch effect removed but biological heterogeneity preserved 
	## sva: estiamte surrogatae variables to remove artifacts

require(sva)


batchQC_num.sv=function(data.matrix, modmatrix){
	data.matrix=as.matrix(data.matrix)
	modmatrix=as.matrix(modmatrix)
    num.sv.be=sva::num.sv(dat=data.matrix,mod=modmatrix, method = "be" , vfilter = NULL, B = 20,seed = 47)
	return(num.sv.be)
}


batchQC_psva=function(data.matrix,batch){
	data.matrix=as.matrix(data.matrix)
	batch=as.factor(batch)
	psva.output=psva(dat=data.matrix, batch=batch)
	return(psva.output)

}

batchQC_sva=function(data.matrix, modmatrix){
    
    n.sv=batchQC_num.sv(data.matrix,modmatrix)
    modmatrix0=model.matrix(~1,data=data.frame(t(data.matrix)))
    sva.object=sva::sva(dat=data.matrix, mod=modmatrix, mod0=modmatrix0, n.sv=n.sv)
    return(sva.object)

}



