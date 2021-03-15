# functions ---------------------------------------------------------------
lgt <- function(x)
  1 / (1 + exp(-x))

rRBM <- function(n, Wr, Wg) {
  X <- matrix(0, n, ncol(Wg), dimnames = list(NULL, colnames(Wg)))
  X[, L[[Q]]] <-
    matrix(rbinom(
      n * length(L[[Q]]),
      size = 1,
      p = matrix(1, n, 1) %x% lgt(matrix(Wr[L[[Q]], 1], 1, length(L[[Q]])))
    ), n, length(L[[Q]]))
  for (i in (Q - 1):2) {
    X[, L[[i]]] <- (X %*% t(Wg))[, L[[i]]]
    X[, L[[i]]] <-
      matrix(rbinom(prod(dim(X[, L[[i]]])), size = 1, c(lgt(X[, L[[i]]]))), dim(X[, L[[i]]])[1], dim(X[, L[[i]]])[2])
  }
  return(X)
}

cdk <- function(X, Wr, Wg) {
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



# cdk_nobias <- function(X, Wr, Wg) {
#   Gr <- Gr2 <- Wr
#   for (i in 3:Q) {
#     X[, L[[i]]] <- (X %*% t(Wr))[, L[[i]]]
#     X[, L[[i]]] <-  matrix(rbinom(prod(dim(X[, L[[i]]])), size = 1, c(lgt(X[, L[[i]]]))), dim(X[, L[[i]]])[1], dim(X[, L[[i]]])[2])
#     # Gr[L[[i]], unlist(L[c( i - 1)])] <-  (t(X) %*% X / n)[L[[i]], unlist(L[c( i - 1)])]
#     
#     Gr[L[[i]], unlist(L[c( i - 1)])] <-  t(X[,L[[i]] ]  ) %*% X[,L[[i-1]] ] / nrow(X)
#     
#   }
#   # Gr[L[[2]], L[[1]]] <- (t(X) %*% X / n)[L[[2]], L[[1]]]
#   X2 <- X
#   for (i in (Q - 1):2) {
#     X2[, L[[i]]] <- (X2 %*% t(Wg))[, L[[i]]]
#     X2[, L[[i]]] <- matrix(rbinom(prod(dim(X2[, L[[i]]])), size = 1, c(lgt(X2[, L[[i]]]))), dim(X2[, L[[i]]])[1], dim(X2[, L[[i]]])[2])
#   }
#   for (i in 3:Q) {
#     X2[, L[[i]]] <-  lgt( (X2 %*% t(Wr))[, L[[i]] ] ) 
#     # layer.lgt <- lgt(X2[, L[[i]]] )
#     # X2[,L[[i]] ]<-layer.lgt
#     if(i<Q) X2[, L[[i]]] <- matrix(rbinom(prod(dim(X2[, L[[i]]])), size = 1, c( X2[, L[[i]]] )), dim(X2[, L[[i]]])[1], dim(X2[, L[[i]]])[2])
#     # X2[, L[[i]]] <-  matrix(rbinom(prod(dim(X2[, L[[i]]])), size = 1, c(lgt(X2[, L[[i]]]))), dim(X2[, L[[i]]])[1], dim(X2[, L[[i]]])[2])
#     
#     Gr2[L[[i]], unlist(L[c( i - 1)])] <-  t(X2[,L[[i]] ]  ) %*% X2[,L[[i-1]] ] / nrow(X2)
#   }
#   # Gr2[L[[2]], L[[1]]] <- (t(X2) %*% X2 / n)[L[[2]], L[[1]]]
#   dG <- Gr - Gr2
#   return(dG)
# }

#Probabilities on the second phase
# cdk_nobias <- function(X, Wr, Wg) {
#   Gr <- Gr2 <- Wr
#   for (i in 3:Q) {
#     X[, L[[i]]] <- (X %*% t(Wr))[, L[[i]]]
#     X[, L[[i]]] <-  matrix(rbinom(prod(dim(X[, L[[i]]])), size = 1, c(lgt(X[, L[[i]]]))), dim(X[, L[[i]]])[1], dim(X[, L[[i]]])[2])
#     # Gr[L[[i]], unlist(L[c( i - 1)])] <-  (t(X) %*% X / n)[L[[i]], unlist(L[c( i - 1)])]
#     
#     Gr[L[[i]], unlist(L[c( i - 1)])] <-  t(X[,L[[i]] ]  ) %*% X[,L[[i-1]] ] / nrow(X)
#     
#   }
#   # Gr[L[[2]], L[[1]]] <- (t(X) %*% X / n)[L[[2]], L[[1]]]
#   X2 <- X
#   for (i in (Q - 1):2) {
#     X2[, L[[i]]] <- (X2 %*% t(Wg))[, L[[i]]]
#     X2[, L[[i]]] <- lgt(X2[, L[[i]]])
#   }
#   for (i in 3:Q) {
#     X2[, L[[i]]] <- (X2 %*% t(Wr))[, L[[i]]]
#     X2[,L[[i]] ] <- lgt(X2[, L[[i]]] )
#     Gr2[L[[i]], unlist(L[c( i - 1)])] <-  t(X2[,L[[i]] ]  ) %*% X2[,L[[i-1]] ] / nrow(X2)
#   }
#   # Gr2[L[[2]], L[[1]]] <- (t(X2) %*% X2 / n)[L[[2]], L[[1]]]
#   dG <- Gr - Gr2
#   return(dG)
# }


genDat<-function(n, p, modelDims){
  
  dat <- matrix(0, n, p)
  q <- length(modelDims)
  dim.total <- 1 + p + sum(modelDims)
  dimlab <- c("bias", paste0("x", 1:p))
  L <- list(dimlab[1], dimlab[-1])
  for (d in 1:q) {
    L[[d + 2]] <- paste0("h", d, "_", 1:modelDims[d])
  }
  dim.v <- unlist(L)
  X <- as.matrix(cbind(1, matrix(0,n,p), matrix(0, n, sum(modelDims))))
  colnames(X) <- dim.v
  W <- matrix(0, dim.total, dim.total,
              dimnames = list(dim.v, dim.v))
  Q <- q + 2
  
  wscale<-2
  for (i in 2:(Q - 1))  W[unlist(L[(i + 1)]), L[[i]]] <-  rnorm(W[unlist(L[(i + 1)]), L[[i]]])*wscale
  
  # for(i in 2:Q) W[L[[i]], L[[1]]] <- sort(rnorm(length(W[L[[i]], L[[1]]]))*wscale)
  for(i in 2:Q) W[L[[i]], L[[1]]] <- seq(-1,1, l=length(W[L[[i]], L[[1]]]))*wscale
  
  
  W[1, 1] <- 1
  Wr <- W
  Wg <- W
  Wg[-1, -1] <- t(Wg[-1, -1])
  Wr[L[[2]] , L[[1]]] <- 0
  #Diags to preser vals
  diag(Wg[L[[Q]], L[[Q]]]) <- 1
  diag(Wr[L[[2]], L[[2]]]) <- 1
  #No bias on first layer
  Wg[L[[Q]], L[[1]]] <- 0
  Wr[L[[1]], L[[2]]] <- 0
  sim.Wg <- Wg
  sim.Wr <- Wr
  # dat <- rRBM(n, sim.Wr, sim.Wg, L)[, 1:ncol(dat) + 1]
  
  X <- matrix(0, n, ncol(Wg), dimnames = list(NULL, colnames(Wg)))
  X[,1]<-1
  X[, L[[Q]]] <-
    matrix(rbinom(
      n * length(L[[Q]]),
      size = 1,
      p = matrix(1, n, 1) %x% lgt(matrix(Wr[L[[Q]], 1], 1, length(L[[Q]])))
    ), n, length(L[[Q]]))
  for (i in (Q - 1):2) {
    X[, L[[i]]] <- (X %*% t(Wg))[, L[[i]]]
    X[, L[[i]]] <-
      matrix(rbinom(prod(dim(X[, L[[i]]])), size = 1, c(lgt(X[, L[[i]]]))), dim(X[, L[[i]]])[1], dim(X[, L[[i]]])[2])
  }
  dat<-X[,L[[2]]]
  # for (i in 3:Q)
  # heatmap.2(sim.Wg[L[[i - 1]], L[[i]]], col = colorRampPalette(c("red", "black", "green")))
  return(list("x"=dat, "Wg"=sim.Wg, "Wr"=sim.Wr))
}


# Pt 2: Simulation --------------------------------------------------------
# Gen data ----------------------------------------------------------------

modelDims<-c(8)
sim<-genDat(20000, 20, modelDims)
dat<-sim$x

sim.Wg<-sim$Wg
sim.Wr<-sim$Wr




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




gc()
closeAllConnections()


lrate = .1
momentum <- 0.8
nIter <- 50000
MSr<-MSg <- 0

Wg.diag <- diag(Wg)
plotFrame<-seq(100,nIter,by=25)
trackVal<-matrix(0, nIter, 4)
for (i in 1:nIter) {
  # cat("\r", i, "\t\t")
  Gr <- cdk(X[sample(1:nrow(X),size=500,replace=F),], Wr, Wg)
  # Gr <- cdk(X[1:100,], Wr, Wg)
  
  Gg <- Gr
  Gg[-1, -1] <- t(Gr)[-1, -1]
  MSr <- MSr * momentum + Gr  
  MSg <- MSg * momentum + Gg
  Wr[Wr != 0] <- Wr[Wr != 0] + MSr[Wr != 0] * lrate
  Wg[Wg != 0] <- Wg[Wg != 0] + MSg[Wg != 0] * lrate
  diag(Wg) <- Wg.diag
  
  trackVal[i,]<- Wg[L[[2]], L[[3]]][1:4] 
  
  
  if(i %in% plotFrame){
    if(i>100)     dev.off()

    par(mfrow=c(1,3))
    plot(trackVal[1:i,1:2], type="o")
    points(trackVal[1:5,1:2], pch=16,cex=2, col="red")
    points(trackVal[(i-5):i,1:2], pch=16,cex=2, col="blue")
    
    plot(trackVal[1:i,2:3], type="o")
    points(trackVal[1:5,2:3], pch=16,cex=2, col="red")
    points(trackVal[(i-5):i,2:3], pch=16,cex=2, col="blue")
    
    plot(trackVal[1:i,3:4], type="o")
    points(trackVal[1:5,3:4], pch=16,cex=2, col="red")
    points(trackVal[(i-5):i,3:4], pch=16,cex=2, col="blue")
    
    
    corMat<-round(cor(sim.Wg[L[[2]],L[[3]]], Wg[L[[2]],L[[3]]]),2)
    corMat.s<-data.frame(corMat)
    corMat.s[abs(corMat)<.5]<-""
    
    corMat2<-round(cor(sim.Wg[L[[3]],L[[4]]], Wg[L[[3]],L[[4]]]),2)
    corMat2.s<-data.frame(corMat2)
    corMat2.s[abs(corMat2)<.5]<-""
    
    # print(corMat.s)
    cat("\014\n")  
    
    cat("\nIteration",i,"\n\n")
    
    print(corMat.s)
    print(corMat2.s)
    
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











visible_state_to_hidden_probabilities <- function(rbm_w, visible_state) {
  1/(1+exp(-rbm_w %*% visible_state))
}

hidden_state_to_visible_probabilities <- function(rbm_w, hidden_state) {
  1/(1+exp(-t(rbm_w) %*% hidden_state))
}

configuration_goodness_gradient <- function(visible_state, hidden_state) {
  hidden_state %*% t(visible_state)/dim(visible_state)[2]
}

sample_bernoulli <- function(mat) {
  dims=dim(mat)
  matrix(rbinom(prod(dims),size=1,prob=c(mat)),dims[1],dims[2])
}

cd1 <- function(rbm_w, visible_data) {
  visible_data = sample_bernoulli(visible_data)
  H0=sample_bernoulli(visible_state_to_hidden_probabilities(rbm_w, visible_data))
  vh0=configuration_goodness_gradient(visible_data, H0)
  V1=sample_bernoulli(hidden_state_to_visible_probabilities(rbm_w, H0))
  H1=visible_state_to_hidden_probabilities(rbm_w, V1)
  vh1=configuration_goodness_gradient(V1, H1)
  vh0-vh1
}

rbm <- function(num_hidden, training_data, learning_rate, n_iterations, mini_batch_size=100, momentum=0.9, quiet=FALSE) {
  #   This trains a model that's defined by a single matrix of weights.
  #   <num_hidden> is the number of hidden units
  #   cd1 is a function that takes parameters <model> and <data> and returns the gradient (or approximate gradient in the case of CD-1) of the function that we're maximizing. Note the contrast with the loss function that we saw in PA3, which we were minimizing. The returned gradient is an array of the same shape as the provided <model> parameter.
  #   This uses mini-batches no weight decay and no early stopping.
  #   This returns the matrix of weights of the trained model.
  n=dim(training_data)[2]
  p=dim(training_data)[1]
  if (n %% mini_batch_size != 0) {
    stop("the number of test cases must be divisable by the mini_batch_size")
  }
  model = (matrix(runif(num_hidden*p),num_hidden,p) * 2 - 1)*.2 
  momentum_speed = matrix(0,num_hidden,p)
  
  start_of_next_mini_batch = 1;
  cat("\n")
  for (iteration_number in 1:n_iterations) {
    if (!quiet) {cat("\r Iter",iteration_number,"\t\t")}
    # mini_batch = training_data[, start_of_next_mini_batch:(start_of_next_mini_batch + mini_batch_size - 1)]
    mini_batch = training_data[,sample(1:ncol(training_data), size=mini_batch_size, replace=FALSE)]
    # start_of_next_mini_batch = (start_of_next_mini_batch + mini_batch_size) %% n
    gradient = cd1(model, mini_batch)
    momentum_speed = momentum * momentum_speed + gradient
    model = model + momentum_speed * learning_rate
  }
  return(model)
}

rbm2<-rbm(num_hidden=modelDims[1], t(dat), learning_rate = .08, n_iterations=1000, mini_batch_size = 100, momentum = .92)


corMat<-round(cor(sim.Wg[L[[2]],L[[3]]], t(rbm2)),2)
corMat.s<-data.frame(abs(corMat)*(abs(corMat)>.7))
corMat.s[corMat.s<.7]<-""
corMat.s


for (i in 3:Q) {
  heatmap.2(sim.Wg[L[[i - 1]], L[[i]]], Rowv=FALSE, col = colorRampPalette(c("red", "black", "green")))
  heatmap.2(t(rbm2), Rowv=FALSE, col = colorRampPalette(c("red", "black", "green")))
}







