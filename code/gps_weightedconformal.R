


source("utils.R")


# -----------------------------------------------------------------------------
# get_gps1
# -----------------------------------------------------------------------------
# Getting propensity scores for spatial interference under conditional ignorability.
# binary treatment, 4 exposure levels 00, 01, 10, 11
# treatment probability ~ multinom or binom

# input
# tr          col vector of individual level binary treatments, numeric
# covar       dataframe/matrix of covariates
# A           spatial network A_ji, in the form of adjacency matrix or distance matrix 
# Atype       type of A, "dist" or "adj", 
#             A does not have to be symmetric. If directed, j in row index
# pstype      "joint" f(x)=P(Z=z,G=g|X=x), 
#             "conditional_g" f(x,z)=P(G=g|Z=z,X=x), 
#             "conditional_z" f(x,g)=P(Z=z|G=g,X=x),
# ps_pred_model  "gbm" or "multinom" or "binomial" 
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
    A0<-A/rowSums(A)
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
    A0<-A/rowSums(A)
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
    } else if (ps_pred_model=="binomial"){
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
    } else if(ps_pred_model=="gbm"){
      T_1 <- gbm::gbm(Y ~ .,distribution="bernoulli",
                      data=data.frame(Y=tr,covar), n.trees = 100)
      T_2 <- gbm::gbm(Y ~ .,distribution="bernoulli",
                      data=data.frame(Y=trnei,covnei), n.trees = 100)               
      ps1 <- predict(T_1, data.frame(covar), type = "response", n.trees = 100)
      ps2 <- predict(T_2, data.frame(covnei), type = "response", n.trees = 100)
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
  } else if (pstype=="conditional_z") {
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

# X, Y
# Xtest     X of testing point. dataframe m rows
# wt_test   wt of testing point. dataframe m rows
# outfun    function for predicting Y
# CQR       whether use quantile regression for Ymodel
# wt        weights of all units 
# alpha     target error rate
# trainid   list/ vector

# output prediction intervals
# cf type="mean", CI two sided


conformalSplit <- function(X, Y, Xtest,wt_test,
                           outfun, CQR,
                           wt,
                           trainid,
                           alpha=0.1) {
  # subset training data and val data
  Xtrain <- X[trainid, ,drop=FALSE]
  Ytrain <- Y[trainid]
  Xcal <- X[-trainid, ,drop=FALSE]
  Ycal <- Y[-trainid]
  
  censoring <- function(x){pmin(pmax(x, 0.05), 20)}
  m <- nrow(Xtest) # the number of testing points
  qt <- c() # vector qt for each test point.
  Yslack <- c()
  
  # need to check weights more closely. why are they behaving so badly?
  # censor first to avoid messing up the mean!
  wt <- pmin(wt,400) 
  wt_test <- pmin(wt_test,400)
  
  wt_cal <- wt[-trainid]
  avg_wt <- mean(wt_cal) # wt_test is a vector
  wt_cal <- censoring(wt_cal/avg_wt)
  totw <- sum(wt_cal)
  wt_cal <- wt_cal / totw
  
  for (i in 1:m) {
    wt_test[i] <- censoring(wt_test[i]/mean(c(wt_cal,wt_test[i])))
    qt[i] <- (1 + wt_test[i] / totw) * (1 - alpha) 
    qt[i] <- pmin(qt[i], 1)
  }
  
  if (CQR==FALSE){
    # mean prediction, default outfun=RF
    #Ymodel <- function(x){do.call(outfun, list(Y = Ytrain, X = Xtrain,Xtest=x))}
    Ymodel <- randomForest::randomForest(x = Xtrain, y = Ytrain)
    Ycal_hat <- predict(Ymodel, newdata = Xcal)
    Yscore <- conformalScore(Ycal, Ycal_hat,quantiles=FALSE) 
    Ytest_hat <- predict(Ymodel, newdata = Xtest) # numeric vector
    
    for (i in 1:m) {Yslack[i]<-weightedConformalCutoff(Yscore, wt_cal, qt[i])} 
    # watch out for the shape of data, Ytest_hat & Yslack both col vec here
    out <-data.frame(Ylo=Ytest_hat,Yup=Ytest_hat)
    out$Ylo<- out$Ylo - Yslack
    out$Yup<- out$Yup + Yslack
    
  } else if (CQR==TRUE){
    # why for CQR Yscore is a vector and not 2 cols?
    if (outfun=="QRF"){
      Ymodel <- grf::quantile_forest(X = Xtrain, Y = Ytrain,quantiles = c(0.05,0.95))
      Ycal_hat <- predict(Ymodel, Xcal,quantiles = c(0.05,0.95))
      Ycal_hat <- Ycal_hat$predictions # matrix col2
      Yscore <- conformalScore(Ycal, Ycal_hat,quantiles=TRUE)
      
      Ytest_hat <- predict(Ymodel, Xtest,quantiles = c(0.05,0.95))
      Ytest_hat <-Ytest_hat$predictions # matrix col2, m rows
      
      for (i in 1:m) {Yslack[i] <- weightedConformalCutoff(Yscore, wt_cal, qt[i]) }
      
      out <-data.frame(Ylo=Ytest_hat[,1],Yup=Ytest_hat[,2])
      out$Ylo<- out$Ylo - Yslack
      out$Yup<- out$Yup + Yslack
      
    } else if (outfun=="QBART"){
      Xtrain <- as.data.frame(Xtrain)
      Ymodel <- bartMachine::bartMachine(X=Xtrain, y=Ytrain, verbose = FALSE)
      
      Ycal_hat <- bartMachine::calc_prediction_intervals(
        Ymodel, new_data = Xcal,pi_conf = 0.95)$interval # matrix col2
      
      Yscore <- conformalScore(Ycal, Ycal_hat,quantiles=TRUE)
      
      Ytest_hat <- bartMachine::calc_prediction_intervals(
        Ymodel, new_data = Xtest,pi_conf = 0.95)$interval # matrix col2, m rows
      
      for (i in 1:m) {Yslack[i] <- weightedConformalCutoff(Yscore, wt_cal, qt[i]) } 
      #Yslack <- weightedquantile(Yscore, prob=0.9, wt_cal,sorted=FALSE) 
      
      out <-data.frame(Ylo=Ytest_hat[,1],Yup=Ytest_hat[,2])
      out$Ylo<- out$Ylo - Yslack
      out$Yup<- out$Yup + Yslack
    }
    
  }
  
  return(out)
}

# -----------------------------------------------------------------------------
# conformalCf_split()
# -----------------------------------------------------------------------------
# conformalCf_split() formats data and params, pass them to conformalSplit() to predict counterfactuals

# to do: take quantiles; take test id list

# X, Y        observed X and Y of all units
# Xtest       data point to do prediction on 
# Ytest       observed Y value for test point
# tltest      observed exposure level of test point 
# wt_test     Propensity of test point receiving observed treatment 
# tl          vector of exposure levels to all units
# A           spatial network A_ji, in the form of adjacency matrix or distance matrix 
# Atype       type of A, "dist" or "adj", 
# ps_pred_model passed to conformalSplit()
# outfun      passed to conformalSplit(), RF/QRF/QBART
# CQR         passed to conformalSplit()
# ps          passed to conformalSplit()
# trainid     optional. interface to other sampling methods. Internal default is random sampling

conformalCf_split <- function(X, Y, 
                              Xtest,Ytest,tltest,wt_test,
                              tl,
                              outfun, CQR,
                              ps,trainid){
  
  inds11 <- which(tl == "t11")
  inds10 <- which(tl == "t10")
  inds01 <- which(tl == "t01")
  inds00 <- which(tl == "t00")
  inds<-list(inds11,inds10,inds01,inds00) # list of index, careful with the ordering!
  
  # parse index
  X0<-Y0<-ps0<-out0<-trainid0<-list()
  n0<-lengths(inds)
  minn<-min(n0)
  trn<-floor(minn*0.6) 
  
  for (i in 1:4){
    if (is.na(trainid)){
      trainid0[[i]]<- sample(n0[[i]], trn) # same size samples
    } else {trainid0[[i]]<-trainid[[i]]}
    X0[[i]] <- X[inds[[i]], ,drop=FALSE]
    Y0[[i]] <- Y[inds[[i]]]
    ps0[[i]]<- ps[inds[[i]], ,drop=FALSE] 
    # X0, Y0 , ps0 has 4 elements, each specific to t11 t10 t01 t00
    # trainid0 elements correspond to each of the X0 Y0 ps0 elements.
    # ps0[[i]] df has 6 (4) cols for 11, 10, 01, 00, (z, g)
  } 
  
  # parse the test points into four groups by tl types.
  # make predictions for each tl type and combine results.
  coln<- c("low11","up11","low10","up10","low01","up01","low00","up00")
  test_grp11<-which(tltest=="t11")# return index within the test points
  if ( length(test_grp11)!=0){
    Xtest_grp11<-Xtest[test_grp11, ,drop=FALSE]
    Ytest_grp11<-Ytest[test_grp11]
    wttest_grp11<-wt_test[test_grp11, ,drop=FALSE]
    
    out1a <- conformalSplit(X0[[1]],Y0[[1]],Xtest_grp11,wttest_grp11[,1]/wttest_grp11[,1],
                            outfun,CQR,wt=ps0[[1]][,1]/ps0[[1]][,1],trainid0[[1]])
    # Tobs is estimated, G could be wrong, therefore Y_Tobs shoudl still be predicted.
    out2a <- conformalSplit(X0[[2]],Y0[[2]],Xtest_grp11,wttest_grp11[,1]/wttest_grp11[,2],
                            outfun,CQR,wt=ps0[[2]][,1]/ps0[[2]][,2],trainid0[[2]])
    
    out3a <- conformalSplit(X0[[3]],Y0[[3]],Xtest_grp11,wttest_grp11[,1]/wttest_grp11[,3],
                            outfun,CQR,wt=ps0[[3]][,1]/ps0[[3]][,3],trainid0[[3]])
    
    out4a <- conformalSplit(X0[[4]],Y0[[4]],Xtest_grp11,wttest_grp11[,1]/wttest_grp11[,4],
                            outfun,CQR,wt=ps0[[4]][,1]/ps0[[4]][,4],trainid0[[4]])
    # out4 cf Y00 is trained on the [[4]] of X0, Y0, ps0, trainid0 where all obs has t00
    
    out11<- cbind.data.frame(out1a,out2a,out3a,out4a)
    colnames(out11)<-coln
    
  } else {
    out11<- data.frame()
    message ("No test points with tl=t11")}
  
  test_grp10<-which(tltest=="t10") 
  if ( length(test_grp10)!=0) {
    Xtest_grp10 <- Xtest[test_grp10, ,drop=FALSE]
    Ytest_grp10 <- Ytest[test_grp10]
    wttest_grp10 <- wt_test[test_grp10, ,drop=FALSE]
    
    out1b <- conformalSplit(X0[[1]],Y0[[1]],Xtest_grp10,wttest_grp10[,2]/wttest_grp10[,1],
                            outfun,CQR,wt=ps0[[1]][,2]/ps0[[1]][,1],trainid0[[1]])
    
    out2b <- conformalSplit(X0[[2]],Y0[[2]],Xtest_grp10,wttest_grp10[,2]/wttest_grp10[,2],
                            outfun,CQR,wt=ps0[[2]][,2]/ps0[[2]][,2],trainid0[[2]])
    
    out3b <- conformalSplit(X0[[3]],Y0[[3]],Xtest_grp10,wttest_grp10[,2]/wttest_grp10[,3],
                            outfun,CQR,wt=ps0[[3]][,2]/ps0[[3]][,3],trainid0[[3]])
    
    out4b <- conformalSplit(X0[[4]],Y0[[4]],Xtest_grp10,wttest_grp10[,2]/wttest_grp10[,4],
                            outfun,CQR,wt=ps0[[4]][,2]/ps0[[4]][,4],trainid0[[4]])
    
    out10<- cbind.data.frame(out1b,out2b,out3b,out4b)
    colnames(out10)<-coln
  } else {
    out10<- data.frame()
    message ("No test points with tl=t10")}
  
  test_grp01 <- which(tltest=="t01") 
  if ( length(test_grp01)!=0) {
    Xtest_grp01 <- Xtest[test_grp01, ,drop=FALSE]
    Ytest_grp01 <- Ytest[test_grp01]
    wttest_grp01 <- wt_test[test_grp01, ,drop=FALSE]
    out1c <- conformalSplit(X0[[1]],Y0[[1]],Xtest_grp01,wttest_grp01[,3]/wttest_grp01[,1],
                            outfun,CQR,wt=ps0[[1]][,3]/ps0[[1]][,1],trainid0[[1]])
    
    out2c <- conformalSplit(X0[[2]],Y0[[2]],Xtest_grp01,wttest_grp01[,3]/wttest_grp01[,2],
                            outfun,CQR,wt=ps0[[2]][,3]/ps0[[2]][,2],trainid0[[2]])
    
    out3c <- conformalSplit(X0[[3]],Y0[[3]],Xtest_grp01,wttest_grp01[,3]/wttest_grp01[,3],
                            outfun,CQR,wt=ps0[[3]][,3]/ps0[[3]][,3],trainid0[[3]])
    
    out4c <- conformalSplit(X0[[4]],Y0[[4]],Xtest_grp01,wttest_grp01[,3]/wttest_grp01[,4],
                            outfun,CQR,wt=ps0[[4]][,3]/ps0[[4]][,4],trainid0[[4]])
    
    out01<- cbind.data.frame(out1c,out2c,out3c,out4c)
    colnames(out01)<-coln
    
  } else {
    out01<- data.frame()
    message ("No test points with tl=t01")}
  
  test_grp00 <- which(tltest=="t00") # checked, correct, numeric vector (integer)
  if ( length(test_grp00)!=0) {
    Xtest_grp00 <- Xtest[test_grp00, ,drop=FALSE]
    Ytest_grp00 <- Ytest[test_grp00]
    wttest_grp00 <- wt_test[test_grp00, ,drop=FALSE]
    out1d <- conformalSplit(X0[[1]],Y0[[1]],Xtest_grp00,wttest_grp00[,4]/wttest_grp00[,1],
                            outfun,CQR,wt=ps0[[1]][,4]/ps0[[1]][,1],trainid0[[1]])
    
    out2d <- conformalSplit(X0[[2]],Y0[[2]],Xtest_grp00,wttest_grp00[,4]/wttest_grp00[,2],
                            outfun,CQR,wt=ps0[[2]][,4]/ps0[[2]][,2],trainid0[[2]])
    
    out3d <- conformalSplit(X0[[3]],Y0[[3]],Xtest_grp00,wttest_grp00[,4]/wttest_grp00[,3],
                            outfun,CQR, wt=ps0[[3]][,4]/ps0[[3]][,3],trainid0[[3]])
    
    out4d <- conformalSplit(X0[[4]],Y0[[4]],Xtest_grp00,wttest_grp00[,4]/wttest_grp00[,4],
                            outfun,CQR,wt=ps0[[4]][,4]/ps0[[4]][,4],trainid0[[4]])
    
    out00<- cbind.data.frame(out1d,out2d,out3d,out4d) 
    colnames(out00)<-coln
    
  } else {
    out00<- data.frame()
    message ("No test points with tl=t00")}
  
  out <- rbind.data.frame(out11,out10,out01,out00)
  out <- na.omit(out)
  # outcome has 8 cols, being the lower and upper bounds of the 4 POs. 
  # outcome has same m of rows as testid length
  return(out)
}

# -----------------------------------------------------------------------------
# conformalCf_split()
# -----------------------------------------------------------------------------
# conformal_ITE1() is the outer func that constructs ITE from Y_obs and counterfactuals

# X, Y
# testid      index of Xtest point in full X. numerical vector.
# tr          individual binary treatments, numeric vector.
# A           spatial network A_ji, in the form of adjacency matrix or distance matrix. 
# Atype       type of A, "dist" or "adj".
# ps_pred_model 
# ps          supply ps if it's known,  diagnostics use only ！
# G           supply G if it's known, diagnostics use only ！
# outfun      passed to conformalSplit()
# CQR         passed to conformalSplit()
# tranid      optional. if supplied, list of length 4. one for each treatment level.
# retrainYmodel If true, Ymodel retrained for each test point. If false, Ymodel fitted only once to make predictions for all points.

# output      ITEs and POs

conformal_ITE1 <- function(X, Y, testid,
                           tr,
                           A, Atype,
                           ps_pred_model="binomial",
                           ps, G,
                           outfun, CQR=FALSE,trainid=NA,
                           retrainYmodel=TRUE
){
  # use psfun="get_gps1"
  PO<-list()
  ITE<-list()
  n <- length(Y)
  m <- length(testid)
  
  # compute joint and marginal propensity scores
  if (is.na(ps)){
    ps<-get_gps1(tr,covar=X,A,Atype,pstype="joint",ps_pred_model) # df n*6/ 4
  } else {ps=ps}
  
  # compute counterfactuals
  if (is.na(G)){
    trnei<-get_trnei1(tr,A,Atype)
  } else {trnei=G}
  
  tl <- paste0("t",tr,trnei) # tls ("t00","t01","t10","t11")
  Xnei<-get_Xnei(covar=X,A=A,Atype=Atype)
  X <- cbind.data.frame(X,Xnei)
  
  if (retrainYmodel==TRUE) {
    for (i in 1:m){
      Ztest <- tr[testid[i]] 
      Gtest <- trnei[testid[i]] 
      tltest <- paste0("t",Ztest,Gtest) 
      wt_test<- ps[testid[i], ,drop=FALSE] # 4/6 cols
      
      #if(tltest=="t11"){wt_test<- ps[testid[i],1]
      #} else if (tltest=="t10") {wt_test<- ps[testid[i],2]
      #} else if (tltest=="t01") {wt_test<- ps[testid[i],3]
      #} else {wt_test<- ps[testid[i],4]}
      
      Xtest <- X[testid[i], ,drop=FALSE]
      Ytest <- Y[testid[i]]
      
      if(is.na(trainid)){
        Cf_params<-list(X,Y,Xtest,Ytest,tltest,wt_test,tl,outfun,CQR,ps=ps,trainid=NA)
      } else {  stop("If trainid not provided, please input NA.")
      }
      # customized trainid under development
      PO_temp<-do.call(conformalCf_split,Cf_params) # get the four POs
      # PO is a df of 8 cols: the lo and up bounds of the 4 POs re 11 10 01 00.
      ITE[[i]]<-data.frame(spill1_lo=PO_temp[[1]] - PO_temp[[4]],
                           spill1_up=PO_temp[[2]] - PO_temp[[3]],
                           spill0_lo=PO_temp[[5]] - PO_temp[[8]],
                           spill0_up=PO_temp[[6]] - PO_temp[[7]],
                           direct1_lo=PO_temp[[1]] - PO_temp[[6]],
                           direct1_up=PO_temp[[2]] - PO_temp[[5]],
                           direct0_lo=PO_temp[[3]] - PO_temp[[8]],
                           direct0_up=PO_temp[[4]] - PO_temp[[7]] )
      PO[[i]]<-PO_temp 
      # to check: is the col order of PO_temp wrong? PO seems to have cols mixed up
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
    PO_temp <-do.call(conformalCf_split,Cf_params) 
    # PO is a df of 8 cols: the lo and up bounds of the 4 POs re 11 10 01 00.
    ITE <- data.frame(spill1_lo=PO_temp[[1]] - PO_temp[[4]],
                      spill1_up=PO_temp[[2]] - PO_temp[[3]],
                      spill0_lo=PO_temp[[5]] - PO_temp[[8]],
                      spill0_up=PO_temp[[6]] - PO_temp[[7]],
                      direct1_lo=PO_temp[[1]] - PO_temp[[6]],
                      direct1_up=PO_temp[[2]] - PO_temp[[5]],
                      direct0_lo=PO_temp[[3]] - PO_temp[[8]],
                      direct0_up=PO_temp[[4]] - PO_temp[[7]] )
    PO[[i]]<-PO_temp
  }
  PO <- do.call(rbind,PO)
  return (list(PO,ITE))
}


