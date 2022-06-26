
# Utils

# -----------------------------------------------------------------------------
# conformal score. for conformalCV main func
# -----------------------------------------------------------------------------
# residual based conformal score
conformalScore <- function(Y, Yhat, quantiles){
  if (quantiles==0) {
    score <- pmax(Yhat - Y, Y - Yhat) # for mean pred
  } else {
    score <- pmax(Yhat[, 1] - Y, Y - Yhat[, 2]) # for CQR
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
# weight conformal cutoff, for conformalCV main func
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
  
  inds <- find_inds(cw, qt)    
  cutoff <- score[inds]
  return(cutoff)
}


# -----------------------------------------------------------------------------
## wrapper for randomForest::randomForest
# -----------------------------------------------------------------------------
RF <- function(Y, X, Xtest){
  # if Y continuous
  fit <- randomForest::randomForest(x = X, y = Y)
  res <- predict(fit, newdata = Xtest)
  return(res)
}


# -----------------------------------------------------------------------------
# quantile random forest. wrapper for grf::quantile_forest 
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
## wrapper for gbm::gbm 
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
  if (Atype=="dist"){
    A0 <- A
    diag(A0) <-1
    A0 <- A^(-1)
    diag(A0) <- 0
    A0 <- A0/rowSums(A0) 
    trnei <- rowSums(matrix(tr,nrow=n,ncol=n,byrow=TRUE)*A0) 
    trnei <- ifelse(trnei>0.5,1,0)
  } else if (Atype=="adj") {
    A0 <- A
    A0 <- A0/rowSums(A0) 
    trnei <- rowSums(matrix(tr,nrow=n,ncol=n,byrow=TRUE)*A0)
    trnei <- ifelse(trnei>0.5,1,0)
  }
  return(trnei)
}