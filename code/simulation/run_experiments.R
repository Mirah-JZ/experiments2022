
source("utils.R")
source("gps_weight_cf.R")


load("datant1_1k5.Rdata")
load("datant1_1k20.Rdata")
load("datant1_5k5.Rdata")
load("datant1_5k20.Rdata")
load("datant2_1k5.Rdata")
load("datant2_1k20.Rdata")
load("datant2_5k5.Rdata")
load("datant2_5k20.Rdata")
load("datant3_1k5.Rdata")
load("datant3_1k20.Rdata")
load("datant3_5k5.Rdata")
load("datant3_5k20.Rdata")
load("datant4_1k5.Rdata")
load("datant4_1k20.Rdata")
load("datant4_5k5.Rdata")
load("datant4_5k20.Rdata")
load("net5k_adj.Rdata")
load("net1k_adj.Rdata")


# summary of coverage and interval lengths

# repeated experiments to calculate coverage
# PO_true is a dataframe of true POs Y11, Y10, Y01, Y00
# return PO coverage rates and lengths of ITE prediction intervals
cal_coverage <- function(X,Y,testid,tr,A,Atype,
                         ps_pred_model,ps,G,
                         outfun,CQR,retrainYmodel,
                         iter=30,PO_true,
                         testID){
  alpha<-c(0.7,0.5,0.3,0.1)
  res<-data.frame(cov11=1:4,cov10=1:4,cov01=1:4,cov00=1:4,
                  len11=1:4,len10=1:4,len01=1:4,len00=1:4)
  for (m in 1:4){
    alpha0 <- alpha[m]
    cov1<-cov2<-cov3<-cov4<-len1<-len2<-len3<-len4<-c()
    for (i in 1:iter){
      temp <- conformal_ITE1(X=X, Y=Y, testid=testid,tr=tr,
                             A=A,Atype=Atype,
                             ps_pred_model=ps_pred_model,
                             ps=ps, G=G,outfun=outfun, CQR=CQR,
                             retrainYmodel=retrainYmodel,
                             alpha=alpha0)
      c1 <- abs(temp[[1]][,2]-temp[[1]][,1]) - pmax(abs(PO_true[,1]-temp[[1]][,1]),
                                                    abs(PO_true[,1]-temp[[1]][,2]))
      c2 <- abs(temp[[1]][,4]-temp[[1]][,3]) - pmax(abs(PO_true[,2]-temp[[1]][,3]),
                                                    abs(PO_true[,2]-temp[[1]][,4]))
      c3 <- abs(temp[[1]][,6]-temp[[1]][,5]) - pmax(abs(PO_true[,3]-temp[[1]][,5]),
                                                    abs(PO_true[,3]-temp[[1]][,6]))
      c4 <- abs(temp[[1]][,8]-temp[[1]][,7]) - pmax(abs(PO_true[,4]-temp[[1]][,7]),
                                                    abs(PO_true[,4]-temp[[1]][,8]))
      cov1[i] <-length(which(c1>0))/length(c1)
      cov2[i] <-length(which(c2>0))/length(c2)
      cov3[i] <-length(which(c3>0))/length(c3)
      cov4[i] <-length(which(c4>0))/length(c4)
      len1[i] <-mean(temp[[2]][,2]-temp[[2]][,1])
      len2[i] <-mean(temp[[2]][,4]-temp[[2]][,3])
      len3[i] <-mean(temp[[2]][,6]-temp[[2]][,5])
      len4[i] <-mean(temp[[2]][,8]-temp[[2]][,7])
    }
    res[m,]<-t(c(mean(cov1),mean(cov2),mean(cov3),mean(cov4),
                 mean(len1),mean(len2),mean(len3),mean(len4)))
  }
  res$mean_cov<-rowMeans(res[,1:4])
  res$mean_len<-rowMeans(res[,6:8])
  res$exp_cov<-c(0.3,0.5,0.7,0.9)
  Y_test<-Y[testid]
  res$exp_meanlen<-c((quantile(Y_test,0.65)-quantile(Y_test,0.35)),
                     (quantile(Y_test,0.75)-quantile(Y_test,0.25)),
                     (quantile(Y_test,0.85)-quantile(Y_test,0.15)),
                     (quantile(Y_test,0.95)-quantile(Y_test,0.05)))
  res$testID<-testID
  
  return(res) # 4*12 dataframe
}

testid_5k <- sample(5000,100) 
testid_1k <- sample(1000,100)
A1<-as.matrix(net1k_adj)
A5<-as.matrix(net5k_adj)


#-------------------------------------------------------------------------------
# test6 linear
#-------------------------------------------------------------------------------
# QBART returns same coverage and 30% shorter intervals, but takes much longer to run
# 0.8 0.8 0.8 1
test6_Ya <- data.frame(Y11=data_nt_1a$Y11[testid_1k],
                       Y10=data_nt_1a$Y10[testid_1k],
                       Y01=data_nt_1a$Y01[testid_1k],
                       Y00=data_nt_1a$Y00[testid_1k])
test6_cov_a <- cal_coverage(X=data_nt_1a$X_nt, Y=data_nt_1a$Y_nt, 
                            testid=testid_1k,
                            tr=data_nt_1a$Z_nt,
                            A=A1,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test6_Ya,
                            testID=61)
Sys.time()
test6_Yb <- data.frame(Y11=data_nt_1b$Y11[testid_1k],
                       Y10=data_nt_1b$Y10[testid_1k],
                       Y01=data_nt_1b$Y01[testid_1k],
                       Y00=data_nt_1b$Y00[testid_1k])
test6_cov_b <- cal_coverage(X=data_nt_1b$X_nt, Y=data_nt_1b$Y_nt, 
                            testid=testid_1k,
                            tr=data_nt_1b$Z_nt,
                            A=A1,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test6_Yb,
                            testID=62)

test6_Yc <- data.frame(Y11=data_nt_1c$Y11[testid_5k],
                       Y10=data_nt_1c$Y10[testid_5k],
                       Y01=data_nt_1c$Y01[testid_5k],
                       Y00=data_nt_1c$Y00[testid_5k])
test6_cov_c <- cal_coverage(X=data_nt_1c$X_nt, Y=data_nt_1c$Y_nt, 
                            testid=testid_5k,
                            tr=data_nt_1c$Z_nt,
                            A=A5,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test6_Yc,
                            testID=63)
Sys.time()
test6_Yd <- data.frame(Y11=data_nt_1d$Y11[testid_5k],
                       Y10=data_nt_1d$Y10[testid_5k],
                       Y01=data_nt_1d$Y01[testid_5k],
                       Y00=data_nt_1d$Y00[testid_5k])
test6_cov_d <- cal_coverage(X=data_nt_1d$X_nt, Y=data_nt_1d$Y_nt, 
                            testid=testid_5k,
                            tr=data_nt_1d$Z_nt,
                            A=A5,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test6_Yd,
                            testID=64)
test6_cov<-rbind.data.frame(test6_cov_a,test6_cov_b,test6_cov_c,test6_cov_d)
save(test6_cov,file="test6_cov_4.Rdata")

#-------------------------------------------------------------------------------
# test7 
#-------------------------------------------------------------------------------
test7_Ya <- data.frame(Y11=data_nt_2a$Y11[testid_1k],
                       Y10=data_nt_2a$Y10[testid_1k],
                       Y01=data_nt_2a$Y01[testid_1k],
                       Y00=data_nt_2a$Y00[testid_1k])
test7_cov_a <- cal_coverage(X=data_nt_2a$X_nt, Y=data_nt_2a$Y_nt, 
                            testid=testid_1k,
                            tr=data_nt_2a$Z_nt,
                            A=A1,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test7_Ya,
                            testID=71)
test7_Yb <- data.frame(Y11=data_nt_2b$Y11[testid_1k],
                       Y10=data_nt_2b$Y10[testid_1k],
                       Y01=data_nt_2b$Y01[testid_1k],
                       Y00=data_nt_2b$Y00[testid_1k])
test7_cov_b <- cal_coverage(X=data_nt_2b$X_nt, Y=data_nt_2b$Y_nt, 
                            testid=testid_1k,
                            tr=data_nt_2b$Z_nt,
                            A=A1,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test7_Yb,
                            testID=72)
test7_Yc <- data.frame(Y11=data_nt_2c$Y11[testid_5k],
                       Y10=data_nt_2c$Y10[testid_5k],
                       Y01=data_nt_2c$Y01[testid_5k],
                       Y00=data_nt_2c$Y00[testid_5k])
test7_cov_c <- cal_coverage(X=data_nt_2c$X_nt, Y=data_nt_2c$Y_nt, 
                            testid=testid_5k,
                            tr=data_nt_2c$Z_nt,
                            A=A5,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test7_Yc,
                            testID=73)
test7_Yd <- data.frame(Y11=data_nt_2d$Y11[testid_5k],
                       Y10=data_nt_2d$Y10[testid_5k],
                       Y01=data_nt_2d$Y01[testid_5k],
                       Y00=data_nt_2d$Y00[testid_5k])
test7_cov_d <- cal_coverage(X=data_nt_2d$X_nt, Y=data_nt_2d$Y_nt, 
                            testid=testid_5k,
                            tr=data_nt_2d$Z_nt,
                            A=A5,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test7_Yd,
                            testID=74)
test7_cov<-rbind.data.frame(test7_cov_a,test7_cov_b,test7_cov_c,test7_cov_d)
save(test7_cov,file="test7_cov_4.Rdata")


#-------------------------------------------------------------------------------
# test8 smoothed X
#-------------------------------------------------------------------------------
test8_Ya <- data.frame(Y11=data_nt_4a$Y11[testid_1k],
                       Y10=data_nt_4a$Y10[testid_1k],
                       Y01=data_nt_4a$Y01[testid_1k],
                       Y00=data_nt_4a$Y00[testid_1k])
test8_cov_a <- cal_coverage(X=data_nt_4a$X_nt, Y=data_nt_4a$Y_nt, 
                            testid=testid_1k,
                            tr=data_nt_4a$Z_nt,
                            A=A1,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test8_Ya,
                            testID=81)
test8_Yb <- data.frame(Y11=data_nt_4b$Y11[testid_1k],
                       Y10=data_nt_4b$Y10[testid_1k],
                       Y01=data_nt_4b$Y01[testid_1k],
                       Y00=data_nt_4b$Y00[testid_1k])
test8_cov_b <- cal_coverage(X=data_nt_4b$X_nt, Y=data_nt_4b$Y_nt, 
                            testid=testid_1k,
                            tr=data_nt_4b$Z_nt,
                            A=A1,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test8_Yb,
                            testID=82)
test8_Yc <- data.frame(Y11=data_nt_4c$Y11[testid_5k],
                       Y10=data_nt_4c$Y10[testid_5k],
                       Y01=data_nt_4c$Y01[testid_5k],
                       Y00=data_nt_4c$Y00[testid_5k])
test8_cov_c <- cal_coverage(X=data_nt_4c$X_nt, Y=data_nt_4c$Y_nt, 
                            testid=testid_5k,
                            tr=data_nt_4c$Z_nt,
                            A=A5,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test8_Yc,
                            testID=83)
test8_Yd <- data.frame(Y11=data_nt_4d$Y11[testid_5k],
                       Y10=data_nt_4d$Y10[testid_5k],
                       Y01=data_nt_4d$Y01[testid_5k],
                       Y00=data_nt_4d$Y00[testid_5k])
test8_cov_d <- cal_coverage(X=data_nt_4d$X_nt, Y=data_nt_4d$Y_nt, 
                            testid=testid_5k,
                            tr=data_nt_4d$Z_nt,
                            A=A5,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test8_Yd,
                            testID=84)
test8_cov<-rbind.data.frame(test8_cov_a,test8_cov_b,test8_cov_c,test8_cov_d)
save(test8_cov,file="test8_cov_4.Rdata")


#-------------------------------------------------------------------------------
# test9 smoothed X, unmeasured confounding
#-------------------------------------------------------------------------------
# remove X1 from test8 X all else the same

test9_cov_a <- cal_coverage(X=data_nt_4a$X_nt[,-1], Y=data_nt_4a$Y_nt, 
                            testid=testid_1k,
                            tr=data_nt_4a$Z_nt,
                            A=A1,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test8_Ya,
                            testID=91)

test9_cov_b <- cal_coverage(X=data_nt_4b$X_nt[,-1], Y=data_nt_4b$Y_nt, 
                            testid=testid_1k,
                            tr=data_nt_4b$Z_nt,
                            A=A1,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test8_Yb,
                            testID=92)

test9_cov_c <- cal_coverage(X=data_nt_4c$X_nt[,-1], Y=data_nt_4c$Y_nt, 
                            testid=testid_5k,
                            tr=data_nt_4c$Z_nt,
                            A=A5,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test8_Yc,
                            testID=93)

test9_cov_d <- cal_coverage(X=data_nt_4d$X_nt[,-1], Y=data_nt_4d$Y_nt, 
                            testid=testid_5k,
                            tr=data_nt_4d$Z_nt,
                            A=A5,Atype="adj",
                            ps_pred_model="multinom",
                            ps=NA, G=NA,
                            outfun="QBART", CQR=TRUE,
                            retrainYmodel=TRUE,
                            iter=30,PO_true=test8_Yd,
                            testID=94)
test9_cov<-rbind.data.frame(test9_cov_a,test9_cov_b,test9_cov_c,test9_cov_d)
save(test9_cov,file="test9_cov_4.Rdata")

#-------------------------------------------------------------------------------
# test10 non-stationary 
#-------------------------------------------------------------------------------
test10_Ya <- data.frame(Y11=data_nt_3a$Y11[testid_1k],
                        Y10=data_nt_3a$Y10[testid_1k],
                        Y01=data_nt_3a$Y01[testid_1k],
                        Y00=data_nt_3a$Y00[testid_1k])
test10_cov_a <- cal_coverage(X=data_nt_3a$X_nt, Y=data_nt_3a$Y_nt, 
                             testid=testid_1k,
                             tr=data_nt_3a$Z_nt,
                             A=A1,Atype="adj",
                             ps_pred_model="multinom",
                             ps=NA, G=NA,
                             outfun="QBART", CQR=TRUE,
                             retrainYmodel=TRUE,
                             iter=30,PO_true=test10_Ya,
                             testID=101)
test10_Yb <- data.frame(Y11=data_nt_3b$Y11[testid_1k],
                        Y10=data_nt_3b$Y10[testid_1k],
                        Y01=data_nt_3b$Y01[testid_1k],
                        Y00=data_nt_3b$Y00[testid_1k])
test10_cov_b <- cal_coverage(X=data_nt_3b$X_nt, Y=data_nt_3b$Y_nt, 
                             testid=testid_1k,
                             tr=data_nt_3b$Z_nt,
                             A=A1,Atype="adj",
                             ps_pred_model="multinom",
                             ps=NA, G=NA,
                             outfun="QBART", CQR=TRUE,
                             retrainYmodel=TRUE,
                             iter=30,PO_true=test10_Yb,
                             testID=102)
test10_Yc <- data.frame(Y11=data_nt_3c$Y11[testid_5k],
                        Y10=data_nt_3c$Y10[testid_5k],
                        Y01=data_nt_3c$Y01[testid_5k],
                        Y00=data_nt_3c$Y00[testid_5k])
test10_cov_c <- cal_coverage(X=data_nt_3c$X_nt, Y=data_nt_3c$Y_nt, 
                             testid=testid_5k,
                             tr=data_nt_3c$Z_nt,
                             A=A5,Atype="adj",
                             ps_pred_model="multinom",
                             ps=NA, G=NA,
                             outfun="QBART", CQR=TRUE,
                             retrainYmodel=TRUE,
                             iter=30,PO_true=test10_Yc,
                             testID=103)
test10_Yd <- data.frame(Y11=data_nt_3d$Y11[testid_5k],
                        Y10=data_nt_3d$Y10[testid_5k],
                        Y01=data_nt_3d$Y01[testid_5k],
                        Y00=data_nt_3d$Y00[testid_5k])
test10_cov_d <- cal_coverage(X=data_nt_3d$X_nt, Y=data_nt_3d$Y_nt, 
                             testid=testid_5k,
                             tr=data_nt_3d$Z_nt,
                             A=A5,Atype="adj",
                             ps_pred_model="multinom",
                             ps=NA, G=NA,
                             outfun="QBART", CQR=TRUE,
                             retrainYmodel=TRUE,
                             iter=30,PO_true=test10_Yd,
                             testID=104)
test10_cov<-rbind.data.frame(test10_cov_a,test10_cov_b,test10_cov_c,test10_cov_d)
save(test10_cov,file="test10_cov_4.Rdata")

#-------------------------------------------------------------------------------
## sensitivity to ps estimation accuracy
#-------------------------------------------------------------------------------
# test 7 1k5
test7_Ya <- data.frame(Y11=data_nt_2a$Y11[testid_1k],
                       Y10=data_nt_2a$Y10[testid_1k],
                       Y01=data_nt_2a$Y01[testid_1k],
                       Y00=data_nt_2a$Y00[testid_1k])
ps7a<-data_nt_2a$ps_nt
ps7a$t11<-ps7a$z*ps7a$g
ps7a$t10<-ps7a$z*(1-ps7a$g)
ps7a$t01<-(1-ps7a$z)*ps7a$g
ps7a$t00<-(1-ps7a$z)*(1-ps7a$g)
ps7a<-ps7a[,-c(1,2)]
nps7a<-ps7a
nps7a$t11<-1;nps7a$t10<-1;nps7a$t01<-1;nps7a$t00<-1
test7_cov_atps <- cal_coverage(X=data_nt_2a$X_nt, Y=data_nt_2a$Y_nt, 
                               testid=testid_1k,
                               tr=data_nt_2a$Z_nt,
                               A=A1,Atype="adj",
                               ps_pred_model="multinom",
                               ps=ps7a, G=NA,
                               outfun="QBART", CQR=TRUE,
                               retrainYmodel=TRUE,
                               iter=10,PO_true=test7_Ya,
                               testID=71)
test7_cov_anps <- cal_coverage(X=data_nt_2a$X_nt, Y=data_nt_2a$Y_nt, 
                               testid=testid_1k,
                               tr=data_nt_2a$Z_nt,
                               A=A1,Atype="adj",
                               ps_pred_model="multinom",
                               ps=nps7a, G=NA,
                               outfun="QBART", CQR=TRUE,
                               retrainYmodel=TRUE,
                               iter=10,PO_true=test7_Ya,
                               testID=71)
test7_cov_aeps <- cal_coverage(X=data_nt_2a$X_nt, Y=data_nt_2a$Y_nt, 
                               testid=testid_1k,
                               tr=data_nt_2a$Z_nt,
                               A=A1,Atype="adj",
                               ps_pred_model="multinom",
                               ps=NA, G=NA,
                               outfun="QBART", CQR=TRUE,
                               retrainYmodel=TRUE,
                               iter=10,PO_true=test7_Ya,
                               testID=71)
test7_cov_atps$ps<-"true"
test7_cov_anps$ps<-"null"
test7_cov_aeps$ps<-"est"
ps_sens<-rbind.data.frame(test7_cov_atps,test7_cov_anps,test7_cov_aeps)

