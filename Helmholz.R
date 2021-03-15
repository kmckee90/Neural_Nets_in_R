

# -------------------------------------------------------------------------
#Helmholtz Machine simulation and training 
#Author: Kevin McKee
#Summary: This generates data from a binary Restricted Boltzmann Machine neural network, which describes the joint distribution of a bit vector using hidden layers.
# The helmholtz machine learns that data generating structure using stochastic gradient descent with a local delta rule by simulating the data.
# i.e., learning happens as a feedback cycle between the recognition model and the generative model.
# The solution should therefore recover the data generating parameters.

# The real-time output shows the correlation between true node weight vectors and currently estimated node weight vectors in layer 1 and 2.
# Layer 2 will not start to match the true values until Layer 1 is mostly correct.

# Some code is left commented to graphically display the learning process of 3 parameters.

# Note: I may code this in RCPP for speed increases of multiple orders of magnitude. Unrepresentatively slow in R.
# -------------------------------------------------------------------------



# functions ---------------------------------------------------------------
lgt<-function(x) 1/(1+exp(-x))

#Generate data top-down (RBM is restricted boltzmann machine)
rRBM <- function(n, Wg) {
  X<-matrix(0, n, ncol(Wg), dimnames=list(NULL,colnames(Wg)) )
  X[,1]<-1
  for(i in Q:2 ){
    X[, L[[i]] ] <- (X %*% t(Wg))[, L[[i]] ]
    X[, L[[i]] ] <- matrix( rbinom(prod(dim(X[, L[[i]] ])), size=1, c(lgt(X[, L[[i]] ])) ),dim(X[, L[[i]] ])[1],dim(X[, L[[i]] ])[2] )
  }
  return(X)
}

#Predict from helmholtz machine: calculate expectations from bottom up then from top down. 
HH.pred <- function(x, Wr, Wg) {
  X<-matrix(0, nrow(x), ncol(Wg), dimnames=list(NULL,colnames(Wg)) )
  X[,1]<-1
  x[,L[[2]]]<-x
  for(i in 3:Q ){
    X[, L[[i]] ] <- lgt(X %*% t(Wr))[, L[[i]] ]
  }
  for(i in (Q-1):2 ){
    X[, L[[i]] ] <- lgt(X %*% t(Wg))[, L[[i]] ]
  }
  return(X)
}


#Set up the dimensions of the simulation, 
dat<-matrix(0, 100000, 100) #Data dimensions
modelDims<-c(30, 10, 5) #Nodes per layer, number of layers given by length of this vector.
BIAS<-T #Include biases in both simulation and fitting?

chunkSize<-30 #Model learns in simple random samples of the data of this size. Small samples tend toward more stochastic learning, more global fit.
lr.s<-.3 #chunkSize/100 seems to work well
maxIter<-900000



# Pt 2: Simulation --------------------------------------------------------
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




# Helmholtz ---------------------------------------------------------------

#Optional annealing i.e., reduction of learning rate
CoolingRate<-0
temp<-(1-.1)*exp(-CoolingRate*(1:maxIter)/maxIter)+.1
plot(temp, type="l")

#Start values
Wr.f<-matrix(F, dim(Wg)[1], dim(Wg)[1], dimnames=dimnames(Wg))
for(i in Q:3) Wr.f[ L[[i]], L[[i-1]] ]<-T
Wg.f<-matrix(F, dim(Wg)[1], dim(Wg)[1], dimnames=dimnames(Wg))
for(i in 2:(Q-1)) Wg.f[ L[[i]], L[[i+1]] ]<-T

if(BIAS){
  Wg.f[unlist(L[2:(Q-1)]),1]<-T
  Wr.f[unlist(L[3:Q]),1]<-T
  Wg.f[-1,1]<-T
  Wr.f[-1,1]<-T
}

Wr[Wr.f]<-rnorm(sum(Wr.f))*.2
Wg[Wg.f]<-rnorm(sum(Wg.f))*.2
diag(Wg[L[[Q]],L[[Q]]])<-0
diag(Wr[L[[2]],L[[2]]])<-0

startPlot<-100
plotFrame<-seq(startPlot,maxIter,by=startPlot)
trackVal<-matrix(NA, maxIter, 4)
  


X<-matrix(0, chunkSize, ncol(Wg), dimnames=list(NULL,colnames(Wg)) )
X[,1]<-1
chunkInd<-1
for(iter in 1:maxIter){

  dat.chunk<-dat[sample(1:nrow(dat), size=chunkSize, replace=FALSE),]
  X[,L[[2]]]<-dat.chunk

  lr<-lr.s*temp[iter]
  for(i in 3:Q){
    X[, L[[i]] ] <- lgt(X %*% t(Wr))[, L[[i]] ]
    X[, L[[i]] ] <- matrix( rbinom(prod(dim(X[, L[[i]] ])), size=1, c(X[, L[[i]] ]) ),dim(X[, L[[i]] ])[1],dim(X[, L[[i]] ])[2] )
  }
  X.p <- lgt(X %*% t(Wg))
  X.p[,1]<-1
  Wg.gr<- t((t(X)%*%(X-X.p))/nrow(X))
  Wg[Wg.f]<- Wg[Wg.f]  +  lr*Wg.gr[Wg.f]

  # SLEEP -------------------------------------------------------------------
  X2<-X
  for(i in Q:2 ){
    X2[, L[[i]] ] <- lgt(X2 %*% t(Wg))[, L[[i]] ]
    X2[, L[[i]] ] <- matrix( rbinom(prod(dim(X2[, L[[i]] ])), size=1, c(X2[, L[[i]] ]) ),dim(X2[, L[[i]] ])[1],dim(X2[, L[[i]] ])[2] )
  }
  X.p <- lgt(X2 %*% t(Wr))
  X.p[,1]<-1
  Wr.gr<- t((t(X2)%*%(X2-X.p))/nrow(X2))
  Wr[Wr.f]<- Wr[Wr.f]  +  lr*Wr.gr[Wr.f]
  


# Plotting ----------------------------------------------------------------
  trackVal[iter,]<- Wg[L[[2]], L[[3]]][1:4] 
  if(iter %in% plotFrame){
    
    #Optional plotting
    # if(iter>(startPlot+1))     dev.off()
    # par(mfrow=c(1,3))
    # plot(trackVal[1:iter,1:2], pch='.')
    # points(trackVal[1:3,1:2], pch=16,cex=2, col="red")
    # points(trackVal[(iter-3):iter,1:2], pch=16,cex=2, col="blue")
    # plot(trackVal[1:iter,2:3], pch='.')
    # points(trackVal[1:3,2:3], pch=16,cex=2, col="red")
    # points(trackVal[(iter-3):iter,2:3], pch=16,cex=2, col="blue")
    # plot(trackVal[1:iter,3:4], pch='.')
    # points(trackVal[1:3,3:4], pch=16,cex=2, col="red")
    # points(trackVal[(iter-3):iter,3:4], pch=16,cex=2, col="blue")
    
    
    corMat<-round(cor(sim.Wg[L[[2]],unlist(L[c(1,3)])], Wg[L[[2]],unlist(L[c(1,3)])]),2)
    corMat.s<-data.frame(corMat)
    corMat.s[abs(corMat)<.7]<-""
    
    cat("\014\n")
    cat("\nIteration",iter,"\n")
    cat("\nLearning rate:",lr,"\n")
    
    cat("\nHidden Layer 1 and 2 weight vector correlations:",lr,"\n")
    

    print(corMat.s)
    
    sim.Wg.0<-sim.Wg[unlist(L[c(2)]),unlist(L[c(1,3)])][,-1]
    est.Wg.0<-Wg[unlist(L[c(2)]),unlist(L[c(1,3)])][,-1]

    dimnames(est.Wg.0)[[2]]<-paste0("est_",dimnames(est.Wg.0)[[2]])

    est.cor<-cor(sim.Wg.0, est.Wg.0)
    est.perm<-ifelse(est.cor>.7,1,est.cor)
    est.perm<-ifelse(est.perm< -.7, -1,est.perm)
    est.perm[abs(est.perm)!=1]<-0
    sim.Wg.1<-(sim.Wg[unlist(L[c(3)]),unlist(L[c(4)])])
    est.Wg.1<-(Wg[unlist(L[c(3)]),unlist(L[c(4)])])
    est.Wg.1.p<-(est.perm)%*%(est.Wg.1)
    
    corMat2<-round(cor(sim.Wg.1, est.Wg.1.p),2)
    corMat2.s<-data.frame(corMat2)
    corMat2.s[abs(corMat2)<.7]<-""
    print(corMat2.s)
    
    
  }
}

sim.Wg.0
round(est.Wg.0,2)

sim.Wg.1
est.Wg.1.p
