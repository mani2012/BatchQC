## library("mlogit")
#library("boot")
#library("pscl")
library("MCMCpack")

##
## Generate simulated count data with batch effects for ngenes
##
seed <- 1000 # For consistent seed across runs
set.seed(seed, kind = NULL, normal.kind = NULL)
nbatch <- 3
ngenes <- 50
prob <- 0.1
#n = 1000
n = 10
mu <- seq(700, length.out=ngenes, by=5)
bmu <- seq(100, length.out=nbatch, by=15000)
bsize <- seq(10, length.out=nbatch, by=1.5)
A.matrix <- matrix(0, nrow=ngenes, ncol=nbatch*n)
for (i in 1:ngenes)  {
  genei <- c()
  for (j in 1:nbatch)  {
    #size = 100
    batchmu <- rnorm(1, mean=bmu[j], sd=1)
    size <- rinvgamma(1, shape=bsize[j], scale=1)
    genei <- c(genei, rnbinom(n,size=size,mu=mu[i]+batchmu))
  }
  A.matrix[i,] <- genei
}

# Heatmap
#cc <- rainbow(ncol(A.matrix))
bc <- rainbow(nbatch)
cc <- rep(bc, each=n)
heatmap(A.matrix, ColSideColors=cc)  

#PCA
pca <- prcomp(t(A.matrix), retx=T, center=T, scale=T) 
pc <- pca$x 
plot(pc[,1], pc[,2],col=cc)

batch <- rep(1:nbatch, each=n)
#Test PCAs
summary(lm(pc[,1]~batch))
summary(lm(pc[,2]~batch))
summary(lm(pc[,3]~batch))
summary(lm(pc[,4]~batch))
summary(lm(pc[,5]~batch))

library(sva)
library(limma)
### apply combat
sample <- 1:nbatch*n
pdata <- data.frame(sample, batch)
modcombat = model.matrix(~1, data=pdata)  
combat_A.matrix = ComBat(dat=A.matrix, batch=batch, mod=modcombat)

# Heatmap
#cc <- rainbow(ncol(combat_A.matrix))
bc <- rainbow(nbatch)
cc <- rep(bc, each=n)
heatmap(combat_A.matrix, ColSideColors=cc)  

#PCA
pca <- prcomp(t(combat_A.matrix), retx=T, center=T, scale=T) 
pc <- pca$x 
plot(pc[,1], pc[,2],col=cc)

batch <- rep(1:nbatch, each=n)
#Test PCAs
summary(lm(pc[,1]~batch))
summary(lm(pc[,2]~batch))
summary(lm(pc[,3]~batch))
summary(lm(pc[,4]~batch))
summary(lm(pc[,5]~batch))

