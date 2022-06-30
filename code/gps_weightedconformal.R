


source("utils.R")


# -----------------------------------------------------------------------------
# get_gps1
# -----------------------------------------------------------------------------
# Getting propensity scores for spatial interference under conditional ignorability.
# Binary treatment, 4 exposure levels 00, 01, 10, 11

# input
# tr          col vector of individual level binary treatments, numeric
# covar       dataframe/matrix of covariates
# A           spatial network A_ji, in the form of adjacency matrix or distance matrix 
# Atype       type of A, "dist" or "adj", 
#             A does not have to be symmetric. If directed, j in row index
# pstype      "joint" f(x)=P(Z=z,G=g|X=x), 
#             "conditional_g" f(x,z)=P(G=g|Z=z,X=x), 
#             "conditional_z" f(x,g)=P(Z=z|G=g,X=x),
# ps_pred_model  "xgboost" or "multinom" or "binomial" 
# output
# ps          if pstype=="joint", outputs joint ps(Z,G) and marginals Z, G
#             if pstype=="conditional_g", outputs marginal ps for G. 
#             if pstype=="conditional_z", outputs marginal ps for Z.

get_gps1 <- function(tr,
                     covar,
                     A,
                     Atype="dist",
                     pstype="joint",
                     ps_pred_model="binomial"){
  # check args insert code here
  
  # format treatment and covariates
  n<-length(tr)
  if (Atype=="dist") {
    # define treatment levels tl based on tr and A, 4 treatment levels
    # to do: these lines need to be put into a function  get_trnei1(tr,A,Atype)
    trnei <-get_trnei1(tr,A,Atype)
    tl <- paste0("t",tr,trnei) # treatment level, multinomial.
    # construct neighborhood covariates
    k <- ncol(covar) 
    covnei <- list()
    for (i in 1:k) {
      covnei[[i]] <- rowSums(matrix(covar[,i],nrow=n,ncol=n,byrow=TRUE)*A0)
    }
    covnei <- as.data.frame(do.call(cbind,covnei))
  } else {
    Atype=="adj"
    # define treatment levels tl
    trnei <-get_trnei1(tr,A,Atype)
    tl <- paste0("t",tr,trnei)
    # construct neighborhood covariates
    k <- ncol(covar)
    covnei <- list()
    for (i in 1:k) {
      covnei[[i]] <- rowSums(matrix(covar[,i],nrow=n,ncol=n,byrow=TRUE)*A0)
    }
    covnei <- as.data.frame(do.call(cbind,covnei))
    netA <- igraph::graph_from_adjacency_matrix(A,mode="undirected",weighted=NULL)
    dgr <- rowSums(A) # the node degrees
    btwnss <- betweenness(netA, weights = NA) # betweeness centrality
    covnei <- cbind.data.frame(covnei,dgr,btwnss)
  }
  covar0 <- cbind.data.frame(covar,covnei)
  # fit treatment model and estimate propensity
  if (pstype=="joint"){
    if (ps_pred_model=="multinom") {
      data0=cbind.data.frame(covar0,tl)
      pred0 <- nnet::multinom(tl~.,data0) 
      ps <- predict(pred0,newdata=data0,type="probs")# probs/ class
      # cols are ordered 00, 01, 10, 11
      ps <-ps[,c(4,3,2,1)]
      colnames(ps)<-c("t11","t10","t01","t00")
    } else {
      # "binomial"
      T_1 <- glm.fit(x=covar,y=tr,family=binomial())
      T_2 <- glm.fit(x=covnei,y=trnei,family=binomial())
      ps1 <-  T_1$fitted.values# p(Z|X_z)
      ps2 <-  T_2$fitted.values # p(G|X_g)
      gps1 <- ps1*ps2 # z=1,g=1
      gps2 <- ps1-gps1 # z=1,g=0 
      gps3 <- (1-ps1)*ps2 # z=0,g=1 
      gps4 <- 1-ps1-gps3 # z=0,g=0 
      ps<-cbind(gps1,gps2,gps3,gps4,ps1,ps2)
      colnames(ps)<-c("t11","t10","t01","t00","z","g")
    } 
  } else if (pstype=="conditional_g"){
    # P(G=g|Z=z,X=x)
    if (ps_pred_model=="xgboost"){
      data0=cbind.data.frame(covar0,tr,trnei)
      pred0 <- xgboost::xgboost(data=data0,label=trnei) 
      ps <- as.matrix(predict(pred0,newdata=data0))
    }else{
      pred0 <- glm.fit(x=c(tr,covar0),y=trnei,family=binomial())
      ps <- pred0$fitted.values
    }    
  } else {
    # "conditional_z" P(Z=z|G=g,X=x)
    if (ps_pred_model=="xgboost"){
      data0=cbind.data.frame(covar0,tr,trnei)
      pred0 <- xgboost::xgboost(data=data0,label=tr) 
      ps <- as.matrix(predict(pred0,newdata=data0))
    }else{
      pred0 <- glm.fit(x=c(trnei,covar0),y=tr,family=binomial())
      ps <- pred0$fitted.values
    }    
  }  
  return(ps)
}


# -----------------------------------------------------------------------------
# conformalSplit()
# -----------------------------------------------------------------------------

# conformalSplit() is the base prediction function, a weighted split cf.

# to do 1: fix function to take a testid list

# X, Y
# Xtest     X of testing point. dataframe 1-m rows
# wt_test   wt of testing point. dataframe 1-m rows
# outfun    function for predicting Y
# CQR       whether use quantile regression for Ymodel
# wt        weights of all units 
# alpha     target error rate
# trainid   list/ vector

# output prediction intervals
# cf type="mean", CI two sided


conformalSplit <- function(X, Y, Xtest,wt_test,
                           outfun, CQR=FALSE,
                           wt,
                           trainid,
                           alpha=0.1) {
  # subset training data and val data
  Xtrain <- X[trainid, ,drop=FALSE]
  Ytrain <- Y[trainid]
  Xcal <- X[-trainid, ,drop=FALSE]
  Ycal <- Y[-trainid]
  
  censoring <- function(x){pmin(pmax(x, 0.05), 20)}
  m <- length(Xtest) # the number of testing points
  qt <- c() # vector qt for each test point.
  for (i in 1:m) {
    wt_cal <- wt[-trainid]
    avg_wt <- mean(c(wt_cal, wt_test[i])) # wt_test is a vector
    wt_cal <- censoring(wt_cal/avg_wt)
    wt_test[i] <- censoring(wt_test[i]/avg_wt)
    totw <- sum(wt_cal)
    #wt_cal <- wt_cal / totw
    qt[i] <- (1 + wt_test[i] / totw) * (1 - alpha) # quantile 1-alpha weighted
    qt[i] <- pmin(qt[i], 1)
  }
  
  if (CQR==FALSE){
    # mean prediction, default outfun=RF
    #Ymodel <- function(x){do.call(outfun, list(Y = Ytrain, X = Xtrain,Xtest=x))}
    Ymodel <- randomForest::randomForest(x = Xtrain, y = Ytrain)
    Ycal_hat <- predict(Ymodel, newdata = Xcal)
    Yscore <- conformalScore(Ycal, Ycal_hat,quantiles=0)
    
    Ytest_hat <- predict(Ymodel, newdata = Xtest) # numeric vector
    
    Yslack <- weightedConformalCutoff(Yscore, wt_cal, qt) 
    Ylo <- Ytest_hat - Yslack
    Yup <- Ytest_hat + Yslack
    out <- data.frame(lower = Ylo, upper = Yup)
  } else{
    # conformalised quantile regression. default outfun= QRF
    Ymodel <- grf::quantile_forest(X = Xtrain, Y = Ytrain,quantiles = c(0.05,0.95))
    Ycal_hat <- predict(Ymodel, Xcal,quantiles = c(0.05,0.95))
    Ycal_hat <- Ycal_hat$predictions
    Yscore <- conformalScore(Ycal, Ycal_hat,quantiles=1)
    
    Ytest_hat <- predict(Ymodel, Xtest,quantiles = c(0.05,0.95))
    Ytest_hat <-Ytest_hat$predictions # matrix
    
    Yslack <- weightedConformalCutoff(Yscore, wt_cal, qt) 
    Ylo <- Ytest_hat[,1] - Yslack
    Yup <- Ytest_hat[,2] + Yslack
    out <- data.frame(lower = Ylo, upper = Yup)
  }
  return(out)
}


# -----------------------------------------------------------------------------
# conformalCf_split()
# -----------------------------------------------------------------------------
# conformalCf_split() formats data and pass them to conformalSplit() 
# It handles the conditioning of cf predictions on the observed tl of test units.

# to do: take quantiles; take test id list

# X, Y        observed X and Y of all units
# Xtest       data point to do prediction on, single point 
# Ytest       observed Y value for test point, single point 
# tltest      observed exposure level of test point, single point 
# wt_test     Propensity of test point receiving observed treatment, single point 
# tl          vector of exposure levels to all units
# A           spatial network A_ji, in the form of adjacency matrix or distance matrix 
# Atype       type of A, "dist" or "adj", 
# ps_pred_model passed to conformalSplit()
# outfun      passed to conformalSplit()
# CQR         passed to conformalSplit()
# ps          passed to conformalSplit()
# trainid     optional. interface to other sampling methods. Internal default is random sampling

conformalCf_split <- function(X, Y, 
                              Xtest,Ytest,tltest,wt_test,
                              tl,
                              outfun, CQR=FALSE,
                              ps,trainid){
  
  inds11 <- which(tl == "t11")
  inds10 <- which(tl == "t10")
  inds01 <- which(tl == "t01")
  inds00 <- which(tl == "t00")
  inds<-list(inds11,inds10,inds01,inds00) # list of index, careful with the ordering
  
  # parse index
  X0<-Y0<-ps0<-out0<-trainid0<-list()
  n0<-lengths(inds)
  minn<-min(n0)
  testn<-floor(minn*0.6)
  
  for (i in 1:4){
    if (is.na(trainid)){
      trainid0[[i]]<- sample(n0[[i]], testn) # same size samples
    } else {trainid0[[i]]<-trainid[[i]]}
    X0[[i]] <- X[inds[[i]], ,drop=FALSE]
    Y0[[i]] <- Y[inds[[i]]]
    ps0[[i]]<- ps[inds[[i]], ,drop=FALSE] 
    # X0, Y0 has 4 elements, each specific to t11 t10 t01 t00
    # trainid0 each element is train id for each of the X0 Y0 elements.
    # ps0 has 6 cols, 11, 10, 01, 00, (z, g)
  } 
  
  
  # parse the test points into four groups by tl types.
  # make predictions for each tl type test points and combine the results.
  test_grp11<-which(tltest=="t11")# return index within the test points
  if ( length(test_grp11)!=0){
    Xtest_grp11<-Xtest[test_grp11, ,drop=FALSE]
    Ytest_grp11<-Ytest[test_grp11]
    wttest_grp11<-wt_test[test_grp11]
    
    out1 <- data.frame(lower=Ytest_grp11,upper=Ytest_grp11)
    
    out2 <- conformalSplit(X0[[2]], Y0[[2]],Xtest_grp11,wttest_grp11,outfun,CQR,
                           wt=ps0[[2]][,1]/ps0[[2]][,2],trainid0[[2]])
    
    out3 <- conformalSplit(X0[[3]], Y0[[3]],Xtest_grp11,wwttest_grp11,outfun,CQR,
                           wt=ps0[[3]][,1]/ps0[[3]][,3],trainid0[[3]])
    
    out4 <- conformalSplit(X0[[4]], Y0[[4]],Xtest_grp11,wttest_grp11,outfun,CQR,
                           wt=ps0[[4]][,1]/ps0[[4]][,4],trainid0[[4]])
    
    out11<- cbind.data.frame(out1,out2,out3,out4)
    
  } else {
    out11<- data.frame()
    message ("No test ponits with tl=11")}
  
  test_grp10<-which(tltest=="t10") 
  if ( length(test_grp10)!=0) {
    Xtest_grp10 <- Xtest[test_grp10, ,drop=FALSE]
    Ytest_grp10 <- Ytest[test_grp10]
    wttest_grp10 <- wt_test[test_grp10]
    
    out1 <- conformalSplit(X0[[1]], Y0[[1]], Xtest_grp10,wttest_grp10,outfun,CQR,
                           wt=ps0[[1]][,2]/ps0[[1]][,1],trainid0[[1]])
    
    out2 <- data.frame(lower=Ytest_grp10,upper=Ytest_grp10)
    
    out3 <- conformalSplit(X0[[3]], Y0[[3]],  Xtest_grp10,wttest_grp10,outfun,CQR,
                           wt=ps0[[3]][,2]/ps0[[3]][,3],trainid0[[3]])
    
    out4 <- conformalSplit(X0[[4]], Y0[[4]],  Xtest_grp10,wttest_grp10,outfun,CQR,
                           wt=ps0[[4]][,2]/ps0[[4]][,4],trainid0[[4]])
    
    out10<- cbind.data.frame(out1,out2,out3,out4)
  } else {
    out10<- data.frame()
    message ("No test ponits with tl=t10")}
  
  test_grp01 <- which(tltest=="t01") 
  if ( length(test_grp01)!=0) {
    Xtest_grp01 <- Xtest[test_grp01, ,drop=FALSE]
    Ytest_grp01 <- Ytest[test_grp01]
    wttest_grp01 <- wt_test[test_grp01]
    out1 <- conformalSplit(X0[[1]], Y0[[1]],Xtest_grp01,wttest_grp01,outfun,CQR,
                           wt=ps0[[1]][,3]/ps0[[1]][,1],trainid0[[1]])
    
    out2 <- conformalSplit(X0[[2]], Y0[[2]], Xtest_grp01,wttest_grp01,outfun,CQR,
                           wt=ps0[[2]][,3]/ps0[[2]][,2],trainid0[[2]])
    
    out3 <- data.frame(lower=Ytest_grp01,upper=Ytest_grp01)
    
    out4 <- conformalSplit(X0[[4]], Y0[[4]], Xtest_grp01,wttest_grp01,outfun,CQR,
                           wt=ps0[[4]][,3]/ps0[[4]][,4],trainid0[[4]])
    
    out01<- cbind.data.frame(out1,out2,out3,out4)
    
  } else {
    out01<- data.frame()
    message ("No test ponits with tl=t01")}
  
  test_grp00 <- which(tltest=="t00") # checked, correct, numeric vector (integer)
  if ( length(test_grp00)!=0) {
    Xtest_grp00 <- Xtest[test_grp00, ,drop=FALSE]
    Ytest_grp00 <- Ytest[test_grp00]
    wttest_grp00 <- wt_test[test_grp00]
    out1 <- conformalSplit(X0[[1]], Y0[[1]], Xtest_grp00,wttest_grp00,outfun,CQR,
                           wt=ps0[[1]][,4]/ps0[[1]][,1],trainid0[[1]])
    
    out2 <- conformalSplit(X0[[2]], Y0[[2]], Xtest_grp00,wttest_grp00,outfun,CQR,
                           wt=ps0[[2]][,4]/ps0[[2]][,2],trainid0[[2]])
    
    out3 <- conformalSplit(X0[[3]], Y0[[3]], Xtest_grp00,wttest_grp00,outfun,CQR,
                           wt=ps0[[3]][,4]/ps0[[3]][,3],trainid0[[3]])
    out4 <- data.frame(lower=Ytest_grp00,upper=Ytest_grp00)
    
    out00<- cbind.data.frame(out1,out2,out3,out4) 
    
  } else {
    out00<- data.frame()
    message ("No test ponits with tl=t00")}
  
  out <- rbind.data.frame(out11,out10,out01,out00)
  # outcome has 8 cols, being the lower and upper bounds of the 4 POs. 
  # outcome has same m of rows as testid length
  return(out)
}


# -----------------------------------------------------------------------------
# conformalCf_split()
# -----------------------------------------------------------------------------

# conformal_ITE1() is the outer func that constructs ITEs
# It merges counterfactual outcomes into POs.

# X, Y
# testid      index of Xtest point in full X. a single index or a numerical vector.
# tr          vector of individual binary treatments, numeric.
# A           spatial network A_ji, in the form of adjacency matrix or distance matrix. 
# Atype       type of A, "dist" or "adj".
# ps_pred_model 
# outfun      passed to conformalSplit()
# CQR         passed to conformalSplit()
# tranid      optional. list of length 4. trainid for each treatment level.
# retrainYmodel If true, Ymodel retrained for each test point. If false, Ymodel fitted only once to make predictions for all points.

# output      ITEs and average ITEs over marginal treatment propensities.

conformal_ITE1 <- function(X, Y, testid,
                           tr,
                           A, Atype,
                           ps_pred_model="binomial",
                           outfun, CQR=FALSE,trainid=NA,
                           retrainYmodel=TRUE
){
  # use psfun="get_gps1"
  
  ITE<-list()
  n <- length(Y)
  m <- length(testid)
  
  # compute joint and marginal propensity scores
  # ps are used here to calculate ITEs, and passed on to conformalSplit() as weights
  ps<-get_gps1(tr,covar=X,A,Atype,pstype="joint",ps_pred_model="binomial") # dataframe n*6
  
  # compute counterfactuals
  trnei<-get_trnei1(tr,A,Atype)
  tl <- paste0("t",tr,trnei) # tls ("t00","t01","t10","t11")
  
  if (retrainYmodel==TRUE) {
    for (i in 1:m){
      Ztest <- tr[testid[i]] 
      Gtest <- trnei[testid[i]] 
      tltest <- paste0("t",Ztest,Gtest) 
      if(tltest=="t11"){wt_test<- ps[testid[i],1]
      } else if (tltest=="t10") {wt_test<- ps[testid[i],2]
      } else if (tltest=="t01") {wt_test<- ps[testid[i],3]
      } else {wt_test<- ps[testid[i],4]
      }
      Xtest <- X[testid[i], ,drop=FALSE]
      Ytest <- Y[testid[i]]
      if(is.na(trainid)){
        Cf_params<-list(X,Y,Xtest,Ytest,tltest,wt_test,tl,outfun,CQR,ps=ps,trainid=NA)
      } else {  stop("If trainid not provided, please input NA.")
      }
      PO<-do.call(conformalCf_split,Cf_params) # get the four potential outcomes
      # PO is a dataframe of 8 cols: the lo and up bounds of the 4 POs re 11 10 01 00.
      ITE[[i]]<-data.frame(spill1_lo=PO[[1]]-PO[[4]],
                           spill1_up=PO[[2]]-PO[[3]],
                           spill0_lo=PO[[5]]-PO[[8]],
                           spill0_up=PO[[6]]-PO[[7]],
                           direct1_lo=PO[[1]]-PO[[6]],
                           direct1_up=PO[[2]]-PO[[5]],
                           direct0_lo=PO[[3]]-PO[[8]],
                           direct0_up=PO[[4]]-PO[[7]] )
    }
    ITE <- do.call("rbind", ITE)
    
  } else if (retrainYmodel==FALSE) {
    Ztest <- tr[testid] # vector
    Gtest <- trnei[testid] # vector
    tltest <- paste0("t",Ztest,Gtest)  # vector
    wt_test <- list()
    for (i in 1:m) {
      if(tltest[i]=="t11"){wt_test[i]<- ps[testid[i],1]
      } else if (tltest[i]=="t10") {wt_test[i]<- ps[testid[i],2]
      } else if (tltest[i]=="t01") {wt_test[i]<- ps[testid[i],3]
      } else {wt_test[i]<- ps[testid[i],4]
      }
    } # list wt_test
    Xtest <- X[testid, ,drop=FALSE]
    Ytest <- Y[testid]
    if(is.na(trainid)){
      Cf_params<-list(X,Y,Xtest,Ytest,tltest,wt_test,tl,outfun,CQR,ps=ps,trainid=NA)
    } else { stop("If trainid not provided, please input NA.")
    }
    PO<-do.call(conformalCf_split,Cf_params) 
    # PO is a dataframe of 8 cols: the lo and up bounds of the 4 POs re 11 10 01 00.
    ITE <- data.frame(spill1_lo=PO[[1]]-PO[[4]],
                      spill1_up=PO[[2]]-PO[[3]],
                      spill0_lo=PO[[5]]-PO[[8]],
                      spill0_up=PO[[6]]-PO[[7]],
                      direct1_lo=PO[[1]]-PO[[6]],
                      direct1_up=PO[[2]]-PO[[5]],
                      direct0_lo=PO[[3]]-PO[[8]],
                      direct0_up=PO[[4]]-PO[[7]] )
  }
  
  return (ITE)
}




