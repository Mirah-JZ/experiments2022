
# Utils

# -----------------------------------------------------------------------------
# conformal score.
# -----------------------------------------------------------------------------
# residual based conformal score
conformalScore <- function(Y, Yhat, quantiles=FALSE){
  if (quantiles==FALSE) {
    score <- pmax(Yhat - Y, Y - Yhat) # for mean pred Y is vector
  } else {
    score <- pmax(Yhat[, 1] - Y, Y - Yhat[, 2]) # for CQR Y interval
  }
  return(score)
}




# -----------------------------------------------------------------------------
# Generate a list of indices for cross-validation. 
# for conformalCf_CV & conformalCV main func
# -----------------------------------------------------------------------------
gen_cv_ids <- function(n, nfolds, offset = 0){
  ids <- sample(n, n)
  quo <- floor(n / nfolds)
  resid <- n - quo * nfolds
  idlist <- lapply(1:nfolds, function(i){
    tmp <- (i - 1) * quo + 1:quo
    if (i <= resid){tmp <- c(tmp, quo * nfolds + i)}
    return(ids[tmp] + offset)}
  )
  return(idlist)
}
# test: gen_cv_ids(50,4)


# -----------------------------------------------------------------------------
# weight conformal cutoff
# -----------------------------------------------------------------------------
weightedConformalCutoff <- function(score, weight, qt){
  
  find_inds <- function(a, b){
    n <- length(a)
    b <- b - 1e-12
    ## n + 1 - rank(-c(a, b), ties.method = "first")[-(1:n)] + rank(-b, ties.method = "first")
    rank(c(a, b), ties.method = "first")[-(1:n)] - rank(b, ties.method = "first") + 1
  }
  
  ord <- order(score)
  weight <- weight[ord]
  score <- score[ord]
  cw <- cumsum(weight)
  
  inds <- find_inds(cw, qt) # find i such that cumsum(w)>=qt
  cutoff <- score[inds]
  return(cutoff)
}


# Weighted quantile function from Tibshirani R package implimentation
# score input is the abs Val set residuals
# prob=1-alpha, 
# w is the weight of val and test points, input is already normalised by totw
weightedquantile <- function(score, prob=0.9, w, sorted=FALSE) {
  if (is.null(w)) w = rep(1,length(score))
  if (!sorted) { o = order(score); score = score[o]; w = w[o] }
  i = which(cumsum(w) >= prob) # if w not normalised use w/sum(w) instead of w
  if (length(i)==0) return(Inf) # Can happen with infinite weights
  else return(score[min(i)])
}


# -----------------------------------------------------------------------------
## randomForest
# -----------------------------------------------------------------------------
# if Y is categorical, res=predict(... , type="prob")
RF <- function(Y, X, Xtest){
  # if Y continuous
  fit <- randomForest::randomForest(x = X, y = Y)
  res <- predict(fit, newdata = Xtest)
  return(res)
}


# -----------------------------------------------------------------------------
# quantile random forest.
# -----------------------------------------------------------------------------
# in main func definitions here, by default quantiles = c(0.05, 0.95)
quantRF <- function(Y, X, Xtest, quantiles, ...){
  fit <- grf::quantile_forest(X, Y, quantiles = quantiles, ...)
  res <- predict(fit, Xtest, quantiles = quantiles)
  if (length(quantiles) == 1){
    res <- as.numeric(res)
  } else {
    res <- as.matrix(res)
  }
  return(res)
}


# -----------------------------------------------------------------------------
## Quantile Boosting
# -----------------------------------------------------------------------------
QBoosting <- function(Y, X, Xtest, quantiles, n.trees = 100){
  if (class(X) != "data.frame"){
    X <- as.data.frame(X)
    Xtest <- as.data.frame(Xtest)
    names(Xtest) <- names(X)
  }
  data <- data.frame(Y = Y, X)
  fit <- gbm::gbm(Y ~ ., distribution = list(name = "quantile", alpha = quantiles[1]), 
                  data = data, n.trees = n.trees)
  res <- predict(fit, Xtest, type = "response", n.trees = n.trees)
  if (length(quantiles) == 2){
    fit2 <- gbm::gbm(Y ~ ., distribution = list(name = "quantile", alpha = quantiles[2]), 
                     data = data, n.trees = n.trees)
    res2 <- predict(fit2, Xtest, type = "response", n.trees = n.trees)
    res <- cbind(res, res2)
  }
  return(res)
}


# -----------------------------------------------------------------------------
## Quantile Bart
# -----------------------------------------------------------------------------
QBART <- function(Y, X, Xtest, quantiles,
                      ndpost = 100){
  if (class(X) != "data.frame"){
    X <- as.data.frame(X)
    Xtest <- as.data.frame(Xtest)
    names(Xtest) <- names(X)
  }
  fit <- bartMachine::bartMachine(X, Y, verbose = FALSE)
  res <- bartMachine::calc_prediction_intervals(
      fit, new_data = Xtest,pi_conf = 0.95)$interval
  
  return(res)
}


# -----------------------------------------------------------------------------
## Bart
# -----------------------------------------------------------------------------
BART <- function(Y, X, Xtest,
                 ndpost = 100, ...){
  if (class(X) != "data.frame"){
    X <- as.data.frame(X)
    Xtest <- as.data.frame(Xtest)
    names(Xtest) <- names(X)
  }
  fit <- bartMachine::bartMachine(X, Y, verbose = FALSE)
  res <- predict(fit, Xtest)
  return(res)
}


# -----------------------------------------------------------------------------
## Boosting
# -----------------------------------------------------------------------------
Boosting <- function(Y, X, Xtest, n.trees = 100){
  # X Xtest should be dataframes
  # suppose Y is binomial (in the use case of {/psfun} for binary treatment )
  if (is.factor(Y)){Y <- as.numeric(Y)-1}
  data <- data.frame(Y = Y, X)
  fit <- gbm::gbm(Y ~ ., distribution = "bernoulli", data = data, n.trees = n.trees, ...)
  res <- predict(fit, Xtest, type = "response", n.trees = n.trees)
  
  return(res)
}

# -----------------------------------------------------------------------------
# spectral clustering for spatial CV folds
# -----------------------------------------------------------------------------
library(RSpectra)
# ref greed::spectral()
spectral <- function(X, K) {
  X <- X + Matrix::t(X)
  ij <- Matrix::which(X > 0, arr.ind = T)
  X[ij] <- 1
  nu <- sum(X) / dim(X)[1]
  D <- (Matrix::rowSums(X) + nu)^-0.5
  Ln <- Matrix::sparseMatrix(ij[, 1], ij[, 2], x = D[ij[, 1]] * D[ij[, 2]])
  V <- RSpectra::eigs(Ln, K)
  Xp <- V$vectors
  Xpn <- Xp / sqrt(rowSums(Xp)^2)
  Xpn[rowSums(Xp) == 0, ] <- 0
  km <- kmeans(Xp, K)
  km$cluster
}

# -----------------------------------------------------------------------------
# get_trnei1 get neighborhood treatment level
# -----------------------------------------------------------------------------
get_trnei1 <- function(tr,A,Atype){
  n<-length(tr)
  if (Atype=="dist"){
    A0<-A
    diag(A0) <-1
    A0 <- A^(-1)
    diag(A0) <- 0
    A0 <- A0/rowSums(A0) # row std inverse dist, can try other dist decay
    trnei <- rowSums(matrix(tr,nrow=n,ncol=n,byrow=TRUE)*A0) 
    trnei <- ifelse(trnei>0.5,1,0) # neighbor treatment {0,1}, 0.5 is arbitrary here
  } else if (Atype=="adj") {
    A0 <- A/rowSums(A) 
    trnei <- rowSums(matrix(tr,nrow=n,ncol=n,byrow=TRUE)*A0)
    trnei <- ifelse(trnei>0.5,1,0)
  }
  return(trnei)
}

# -----------------------------------------------------------------------------
# get_Xnei get neighborhood mean X 
# -----------------------------------------------------------------------------
get_Xnei <- function(covar,A,Atype){
  n<-nrow(covar)
  k <- ncol(covar)
  covnei <- list()
  if(Atype=="dist") {
    A0<-A
    diag(A0) <-1
    A0 <- A^(-1)
    diag(A0) <- 0
    A0 <- A0/rowSums(A0) # row std inverse dist
    for (i in 1:k) {
      covnei[[i]] <- rowSums(matrix(covar[,i],nrow=n,ncol=n,byrow=TRUE)*A0)
    }
    covnei <- as.data.frame(do.call(cbind,covnei))
    
  } else if (Atype=="adj") {
    A0 <- A/rowSums(A)
    for (i in 1:k) {
      covnei[[i]] <- rowSums(matrix(covar[,i],nrow=n,ncol=n,byrow=TRUE)*A0)
    }
    covnei <- as.data.frame(do.call(cbind,covnei))
    netA <- igraph::graph_from_adjacency_matrix(A,mode="undirected",weighted=NULL)
    dgr <- rowSums(A) # the node degrees
    btwnss <- betweenness(netA, weights = NA) # betweeness centrality
    covnei <- cbind.data.frame(covnei,dgr,btwnss)
  }
  return(covnei)
}
