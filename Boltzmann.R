
# -------------------------------------------------------------------------
# Code for Restricted B Machine (RBM)
# Author: Kevin McKee
# Code is abstracted to allow any number of layers and nodes per layer.
# Uses simple delta rule to learn.
# Plotting functions are available but commented out.
# Output is currently in the form of a correlation matrix between the true node structures
# and the estimated node structures.
# -------------------------------------------------------------------------




library(gplots)

# functions ---------------------------------------------------------------
lgt <- function(x)
  1 / (1 + exp(-x))


rRBM <- function(n, Wg) {
  X<-matrix(0, n, ncol(Wg), dimnames=list(NULL,colnames(Wg)) )
  X[,1]<-1
  for(i in Q:2 ){
    X[, L[[i]] ] <- (X %*% t(Wg))[, L[[i]] ]
    X[, L[[i]] ] <- matrix( rbinom(prod(dim(X[, L[[i]] ])), size=1, c(lgt(X[, L[[i]] ])) ),dim(X[, L[[i]] ])[1],dim(X[, L[[i]] ])[2] )
  }
  return(X)
}


calcGrad <- function(X, Wr, Wg) {
  Gr <- Gr2 <- Wr
  for (i in 3:Q) {
    X[, L[[i]]] <- (X %*% t(Wr))[, L[[i]]]
    X[, L[[i]]] <-  matrix(rbinom(prod(dim(X[, L[[i]]])), size = 1, c(lgt(X[, L[[i]]]))), dim(X[, L[[i]]])[1], dim(X[, L[[i]]])[2])
    Gr[L[[i]], unlist(L[c(1, i - 1)])] <- (t(X) %*% X / nrow(X))[L[[i]], unlist(L[c(1, i - 1)])]
  }
  Gr[L[[2]], L[[1]]] <- (t(X) %*% X / nrow(X))[L[[2]], L[[1]]]
  X2 <- X
  for (i in (Q - 1):2) {
    X2[, L[[i]]] <- (X2 %*% t(Wg))[, L[[i]]]
    X2[, L[[i]]] <-  matrix(rbinom(prod(dim(X2[, L[[i]]])), size = 1, c(lgt(X2[, L[[i]]]))), dim(X2[, L[[i]]])[1], dim(X2[, L[[i]]])[2])
  }
  for (i in 3:Q) {
    X2[, L[[i]]] <-  lgt( (X2 %*% t(Wr))[, L[[i]] ] ) 
    if(i<Q) X2[, L[[i]]] <- matrix(rbinom(prod(dim(X2[, L[[i]]])), size = 1, c( X2[, L[[i]]] )), dim(X2[, L[[i]]])[1], dim(X2[, L[[i]]])[2])
    Gr2[L[[i]], unlist(L[c(1, i - 1)])] <- (t(X2) %*% X2 / nrow(X))[L[[i]], unlist(L[c(1, i - 1)])]
  }
  Gr2[L[[2]], L[[1]]] <- (t(X2) %*% X2 / nrow(X) )[L[[2]], L[[1]]]
  dG <- Gr - Gr2
  return(dG)
}



# Pt 2: Simulation --------------------------------------------------------
# Gen data ----------------------------------------------------------------

modelDims<-c(8)
BIAS<-TRUE

# Gen data ----------------------------------------------------------------
p<-ncol(dat)
n<-nrow(dat)
q<-length(modelDims)
dim.total<-1+p+sum(modelDims)
dimlab<-c("bias",paste0("x",1:p))
L<-list(dimlab[1], dimlab[-1])
for(d in 1:q) {
  L[[d+2]]<-paste0("h",d,"_",1:modelDims[d])
}
dim.v<-unlist(L)
X<-as.matrix(cbind(1, dat, matrix(0,n,sum(modelDims))))
colnames(X)<-dim.v
W<-matrix(0, dim.total, dim.total,dimnames=list(dim.v,dim.v))
Q<-q+2
for(i in 2:(Q-1) ) W[ unlist( L[ (i+1)]), L[[i]] ] <- rnorm( W[ unlist( L[ (i+1)]), L[[i]] ])*2
W[   unlist( L[ 2:Q]) , L[[1]] ] <- BIAS*rnorm(length(W[   unlist( L[ 2:Q]) , L[[1]] ]))
W[1,1]<-1
Wr<-W
Wg<-W
Wg[-1,-1]<-t(Wg[-1,-1])
Wr[ L[[2]] , L[[1]] ] <- 0
#Diags to preser vals
diag(Wg[L[[Q]],L[[Q]]])<-1
diag(Wr[L[[2]],L[[2]]])<-1
#No bias on first layer
Wg[L[[Q]],L[[1]]]<-0
Wr[L[[1]],L[[2]]]<-0
sim.Wg<-Wg
sim.Wr<-Wr

X.gen<-rRBM(n, sim.Wg)
dat<-X.gen[,L[[2]]]






# Fit to simulated data ---------------------------------------------------
p <- ncol(dat)
n <- nrow(dat)
q <- length(modelDims)
Q <- q + 2
dim.total <- 1 + p + sum(modelDims)
dimlab <- c("bias", paste0("x", 1:p))
L <- list(dimlab[1], dimlab[-1])
for (d in 1:q) {
  L[[d + 2]] <- paste0("h", d, "_", 1:modelDims[d])
}
dim.v <- unlist(L)
X <- as.matrix(cbind(1, dat, matrix(0, n, sum(modelDims))))
colnames(X) <- dim.v



W <- matrix(0, dim.total, dim.total,
            dimnames = list(dim.v, dim.v))
for (i in 2:(Q - 1))  W[unlist(L[(i + 1)]), L[[i]]] <-runif(W[unlist(L[(i + 1)]), L[[i]]])*2-1
W[unlist(L[2:Q]) , L[[1]]] <-  (rnorm(length(W[unlist(L[2:Q]) , L[[1]]])))
W[1, 1] <- 1
Wr<- Wg <- W
Wg[-1, -1] <- t(Wg[-1, -1])
Wr[L[[2]] , L[[1]]] <- 0
#Diags to preser vals
diag(Wg[L[[Q]], L[[Q]]]) <- 1
diag(Wr[L[[2]], L[[2]]]) <- 1
#No bias on first layer
Wg[L[[Q]], L[[1]]] <- 0
Wr[L[[1]], L[[2]]] <- 0



lrate = .1
momentum <- 0.8
nIter <- 50000
MSr<-MSg <- 0

Wg.diag <- diag(Wg)
plotFrame<-seq(100,nIter,by=25)
trackVal<-matrix(0, nIter, 4)
for (i in 1:nIter) {
  # cat("\r", i, "\t\t")
  Gr <- calcGrad(X[sample(1:nrow(X),size=500,replace=F),], Wr, Wg)

  
  Gg <- Gr
  Gg[-1, -1] <- t(Gr)[-1, -1]
  MSr <- MSr * momentum + Gr  
  MSg <- MSg * momentum + Gg
  Wr[Wr != 0] <- Wr[Wr != 0] + MSr[Wr != 0] * lrate
  Wg[Wg != 0] <- Wg[Wg != 0] + MSg[Wg != 0] * lrate
  diag(Wg) <- Wg.diag
  
  trackVal[i,]<- Wg[L[[2]], L[[3]]][1:4] 
  
  # 
  if(i %in% plotFrame){
    
  # if(i>200)    dev.off()
    # par(mfrow=c(1,3))
    # plot(trackVal[1:i,1:2], type="o")
    # points(trackVal[1:5,1:2], pch=16,cex=2, col="red")
    # points(trackVal[(i-5):i,1:2], pch=16,cex=2, col="blue")
    # 
    # plot(trackVal[1:i,2:3], type="o")
    # points(trackVal[1:5,2:3], pch=16,cex=2, col="red")
    # points(trackVal[(i-5):i,2:3], pch=16,cex=2, col="blue")
    # 
    # plot(trackVal[1:i,3:4], type="o")
    # points(trackVal[1:5,3:4], pch=16,cex=2, col="red")
    # points(trackVal[(i-5):i,3:4], pch=16,cex=2, col="blue")


    corMat<-round(cor(sim.Wg[L[[2]],L[[3]]], Wg[L[[2]],L[[3]]]),2)
    corMat.s<-data.frame(corMat)
    corMat.s[abs(corMat)<.5]<-""

    # corMat2<-round(cor(sim.Wg[L[[3]],L[[4]]], Wg[L[[3]],L[[4]]]),2)
    # corMat2.s<-data.frame(corMat2)
    # corMat2.s[abs(corMat2)<.5]<-""
    cat("\014\n")
    
    print(corMat.s)
    # print(corMat2.s)
    cat("\nIteration",i,"\n\n")

  }
}





require(gplots)
#Results
cor(sim.Wg[ L[[2]], 1],  Wg[ L[[2]] ,1])
cor(sim.Wr[ L[[3]], 1],  Wr[ L[[3]] ,1])

heatmap.2(cbind(sim.Wg[-1,1], Wg[-1,1], sim.Wr[-1,1], Wr[-1,1]),Rowv=F, Colv=F, col=colorRampPalette(c("red","black","blue")), main="Bias")
for (i in 3:3) {
  heatmap.2(sim.Wg[L[[i - 1]], L[[i]]], Rowv=FALSE, col = colorRampPalette(c("red", "black", "green")))
  heatmap.2(Wg[L[[i - 1]], L[[i]]], Rowv=FALSE, col = colorRampPalette(c("red", "black", "green")))
}







