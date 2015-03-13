### simulate data
nbatch <- 3
ncond <- 2
npercond <- 10
data.matrix <- rnaseq_sim(ngenes=50, nbatch=nbatch, ncond=ncond, npercond=npercond, 
                          ggstep=5, bbstep=15000, ccstep=10000, bvarstep=2, seed=1234)

### apply BatchQC
batch <- rep(1:nbatch, each=ncond*npercond)
condition <- rep(rep(1:ncond, each=npercond), nbatch)
nsample <- nbatch*ncond*npercond
sample <- 1:nsample
pdata <- data.frame(sample, batch, condition)
modmatrix = model.matrix(~as.factor(condition), data=pdata)  
pca <- batchQC(data.matrix, batch=batch, mod=modmatrix)

### apply combat
combat_data.matrix = ComBat(dat=data.matrix, batch=batch, mod=modmatrix)

### Rerun the BatchQC pipeline on the batch adjusted data
pca <- batchQC(combat_data.matrix, batch=batch, mod=modmatrix)
