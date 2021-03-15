

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


# Real data ---------------------------------------------------------------
setwd("D:/VT/Code/NN/rbm")
dat<-read.csv("testTextDat.csv", header = FALSE)[,1]
# table(dat)
# as.matrix(sort(table(dat)))

p<-ceiling(log2(length(unique(dat))))
bitMap<-expand.grid(rep( list(0:1), p))[1:length(unique(dat)),]
rownames(bitMap)<-sort(unique(dat))
dat.bit<-matrix(NA, length(dat),p)
for(i in 1:length(dat)) dat.bit[i,1:10]<-unlist(bitMap[dat[i],])


dat<-dat.bit












# functions ---------------------------------------------------------------
lgt<-function(x) 1/(1+exp(-x))

#Generate data top-down (RBM is restricted boltzmann machine)
rRHH <- function(n, Wg) {

  X<-matrix(0, n, dim.total, dimnames=list(NULL, dim.v))
  X[,1]<-1
  for(t in 1:n){
    if(t>1) for(l in 1:min(t-1,lags) ) X[t,unlist(lapply(L[-1],"[",l+1))]<-X[t-l,unlist(lapply(L[-1],"[",1))]
    x<-as.matrix(X[t,])
    for(h in Q:2){
      x[L[[h]][[1]],]<-lgt(Wg%*%x)[L[[h]][[1]],]
      x[L[[h]][[1]],]<-rbinom(length(x[L[[h]][[1]]]), 1, x[L[[h]][[1]],] )
    }
    X[t,unlist(lapply(L,"[",1))]<-t(x[unlist(lapply(L,"[",1)),])
    
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
printData<-function(x){
  dat.char<-matrix(NA, nrow(x), 1)
  k1<-apply(x,1, paste0, collapse="")
  k2<-apply(bitMap,1, paste0, collapse="")
  for(i in 1:length(k1)) {
    dat.char[i]<-c(names(k2)[k2==k1[i]], NA)[1]
  }
  return(paste0(dat.char, collapse=" "))
} 


# dat<-matrix(0, 100, 10) #Data dimensions

#Set up the dimensions of the simulation, 
n<-nrow(dat)
p<-ncol(dat)
modelDims<-rep(10, 10) #Nodes per layer, number of layers given by length of this vector.
lags<-10
chunkSize<-20 #Model learns in simple random samples of the data of this size. Small samples tend toward more stochastic learning, more global fit.
lr<-.1
# lr.s<-.3 #chunkSize/100 seems to work well
maxIter<-900000



# Pt 2: Simulation --------------------------------------------------------
# Gen data ----------------------------------------------------------------

# p<-ncol(dat)
# n<-nrow(dat)
q<-length(modelDims)
Q<-q+2

# dim.total<-1+(lags+1)*p+(lags+1)*sum(modelDims)
dimlab<-paste0("x",1:p)
# dimlab<-paste0(rep(dimlab,lags), "_l",rep(0:lags, each=p))

lagList<-list()
for(l in 0:lags){
  lagList<-c(lagList, list(paste0(dimlab,"_l",l)) )
}
L<-list("bias",  lagList)
for(d in 1:q) {
  lagList<-list()
  for(l in 0:lags){  
    hlab<-paste0("h",d,"_n",1:modelDims[d])
    lagList<-c(lagList, list(paste0(hlab,"_l",l )))
  }
  L<-c(L, list(lagList))
}

# dim.v<-unlist(L)
dim.v<-L[[1]]
for(l in 0:lags+1 )for(h in 2:Q) dim.v<-c(dim.v, L[[h]][[l]])


dim.total<-length(dim.v)
X<-matrix(0, n, dim.total, dimnames=list(NULL, dim.v))
X[,1]<-1

t0dims<-unlist(lapply(L[-1],"[",1)) #dim.v[1:(2+sum(modelDims))+1]

W<-matrix(0, dim.total, dim.total,dimnames=list(dim.v,dim.v))


W.f<-W==1
W[1,1]<-1

W.f[t0dims,1]<-TRUE

for(i in 3:Q){
  W.f[L[[i-1]][[1]],L[[i]][[1]]]<-T
}

for(i in 2:Q){
  for(l in 1:(lags))
    W.f[ L[[i]][[1]], L[[i]][[l+1]] ]<-T
    try(W.f[ L[[i]][[1]], L[[i+1]][[l]] ] <- T, silent=T)
}

W.f[t0dims,setdiff(unlist(L), c(t0dims,"bias") )]<-T

W[W.f]<-rnorm(sum(W.f))
Wg<-W

Wr<-Wg
Wr[t0dims,t0dims]<-t(Wg[t0dims,t0dims])
# Wg[]
round(Wg,2);round(Wr,2)
sim.Wg<-Wg

Wg.f<-W.f
Wr.f<-W.f;Wr.f[t0dims,t0dims]<-t(Wg.f[t0dims,t0dims])

  
# X.gen<-rRHH(n, sim.Wg)
# dat<-X.gen[,L[[2]]]
# dat<-X.gen[,unlist(L[[2]][[1]])]




# Helmholtz ---------------------------------------------------------------
#Start values
Wg[Wg.f]<-rnorm(sum(Wg.f))
Wr[Wr.f]<-rnorm(sum(Wr.f))


Wr[Wr==1]<-0
Wg[Wg==1]<-0
# Wr[1,1]<-1
# Wg[1,1]<-1



nIter<- 500000
for(i in 1:nIter){
  # X.gen<-rRHH(n+1, sim.Wg)[-1,]
  # dat<-X.gen[,unlist(L[[2]][[1]])]
  dat.chunk.seg<-sample(1:(nrow(dat)-chunkSize), size=1 )
  dat.chunk<-dat[dat.chunk.seg:(dat.chunk.seg+chunkSize-1),]
  
  n<-chunkSize
  
  X<-matrix(0, n, dim.total, dimnames=list(NULL, dim.v))
  X[,1]<-1
  X[,unlist(L[[2]][[1]]) ]<-dat.chunk
  
  
#WAKE
  for(t in 1:n){
    if(t>1) for(l in 1:min(t-1,lags) ) X[t,unlist(lapply(L[-1],"[",l+1))]<-X[t-l,unlist(lapply(L[-1],"[",1))]
    x<-as.matrix(X[t,])
    for(h in 3:Q){
      x[L[[h]][[1]],]<-lgt(Wr%*%x)[L[[h]][[1]],]
      x[L[[h]][[1]],]<-rbinom(length(x[L[[h]][[1]]]), 1, x[L[[h]][[1]],] )
    }
    X[t,unlist(lapply(L,"[",1))]<-t(x[unlist(lapply(L,"[",1)),])
  }
  
   X.p<-lgt(X%*%t(Wg))
   X.p[,1]<-1
   Wg.gr<-t((t(X)%*%(X-X.p))/nrow(X))
   Wg[Wg.f]<- (Wg + lr*Wg.gr)[Wg.f]
  
    
#SLEEP
   X<-matrix(0, n, dim.total, dimnames=list(NULL, dim.v))
   X[,1]<-1
   for(t in 1:n){
     if(t>1) for(l in 1:min(t-1,lags) ) X[t,unlist(lapply(L[-1],"[",l+1))]<-X[t-l,unlist(lapply(L[-1],"[",1))]
     x<-as.matrix(X[t,])
     for(h in Q:2){
       x[L[[h]][[1]],]<-lgt(Wg%*%x)[L[[h]][[1]],]
       x[L[[h]][[1]],]<-rbinom(length(x[L[[h]][[1]]]), 1, x[L[[h]][[1]],] )
     }
     X[t,unlist(lapply(L,"[",1))]<-t(x[unlist(lapply(L,"[",1)),])
     
   }
  
    
   X.p<-lgt(X%*%t(Wr))
   X.p[,1]<-1
   Wr.gr<-t((t(X)%*%(X-X.p))/nrow(X))
   Wr[Wr.f]<-(Wr + lr*Wr.gr)[Wr.f]
   
   if(i%%5==0) cat("\r",i,"\t\t")
   
   #Output some shit
   if(i%%50==0){
     cat("\nData:\t")
     print(printData(dat.chunk))
     cat("\nSleep output:\t")
     print(printData(X[,L[[2]][[1]]]))
     cat("\n\nSimulation:\t")
     print(printData(rRHH(20, Wg)[,L[[2]][[1]]]))
     cat("\n\n")
   }
}


 


printData(simDat)

# Generate some shit ------------------------------------------------------
simDat<-rRHH(20, Wg)[,L[[2]][[1]]]

dat.char<-matrix(NA, nrow(simDat), 1)

k1<-apply(simDat,1, paste0, collapse="")
k2<-apply(bitMap,1, paste0, collapse="")
for(i in 1:length(k1)) {
  dat.char[i]<-c(names(k2)[k2==k1[i]], NA)[1]
}
print(paste0(dat.char, collapse=" "))
}