
# Utils

# -----------------------------------------------------------------------------
# conformal score. for conformalCV main func
# -----------------------------------------------------------------------------
# residual based conformal score
conformalScore <- function(Y, Yhat){
  score <- pmax(Yhat[, 1] - Y, Y - Yhat[, 2])
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
RF <- function(Y, X, Xtest, ...){
  # depending on the class of Y being continuous or binomial
  # if continuous
  fit <- randomForest::randomForest(x = X, y = Y, ...)
  res <- predict(fit, newdata = Xtest)
  res <- as.numeric(res)
  # can be extended to other cases ...
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
# spectral clustering for spatal CV folds
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