
library(gstat)
library(fields)
library(sp)
library(igraph)
library(tidyverse)
library(ggplot2)



# Generating non-spatial data

## dataset_1 dataset_2
# non-spatial data, linear treatment effect.
gen_ns <- function(n,d,linear=TRUE){
  set.seed(123)
  # X_ns <- matrix(rnorm(n * d), nrow = n) # norm
  X_ns <- matrix(runif(n * d,-3,3), nrow = n, ncol = d) # unif [-3,3]
  beta_ns <- rep(1, d)
  if (linear==TRUE){
    Y1_ns <- X_ns %*% beta_ns + rnorm(n)
  } else {
    genY_ns_2 <- function(X){
      (10/(1+exp(-12*(X[, 1]-0.5))))*(2/(1+exp(-12*(X[, 2]-0.5))))+
        X_ns %*% beta_ns-10+rnorm(n)
    }
    Y1_ns <- genY_ns_2(X_ns)
  }
  
  Y0_ns <- rnorm(n)
  # ps_ns_1 <- (1 + pbeta(X_ns_1[, 1], 2, 4)) / 4 # [0.25,0.5] 
  # T_ns_1 <- as.numeric(runif(n) < ps_ns_1) # 
  ps_ns <- pnorm(X_ns[, 1]) 
  T_ns <- sapply(ps_ns,function(x){rbinom(1,1,prob=x)})
  Y_ns <- ifelse(T_ns == 1, Y1_ns, Y0_ns)
  return(list(X_ns=X_ns,
              Y1_ns=Y1_ns,Y0_ns=Y0_ns,Y_ns=Y_ns,
              ps_ns=ps_ns,T_ns=T_ns)) 
}

data_ns_1a<-gen_ns(1000,5,linear=TRUE)
data_ns_1b<-gen_ns(1000,20,linear=TRUE)  
data_ns_1c<-gen_ns(5000,5,linear=TRUE)  
data_ns_1d<-gen_ns(5000,20,linear=TRUE)  

hist(data_ns_1d$X_ns[,1])
hist(data_ns_1d$ps_ns)
hist(data_ns_1d$Y1_ns)
sd(data_ns_1d$Y1_ns) # 7.83
plot(x=data_ns_1d$X_ns[,1],y=data_ns_1d$Y1_ns)

save(data_ns_1a,file="datans1_1k5.Rdata")
save(data_ns_1b,file="datans1_1k20.Rdata")
save(data_ns_1c,file="datans1_5k5.Rdata")
save(data_ns_1d,file="datans1_5k20.Rdata")

data_ns_2a<-gen_ns(1000,5,linear=FALSE)
data_ns_2b<-gen_ns(1000,20,linear=FALSE)  
data_ns_2c<-gen_ns(5000,5,linear=FALSE)  
data_ns_2d<-gen_ns(5000,20,linear=FALSE)  

hist(data_ns_2a$X_ns[,1])
hist(data_ns_2d$ps_ns)
hist(data_ns_2d$Y1_ns)
sd(data_ns_2d$Y1_ns) # 11.8 more or less on par with linear data
plot(x=data_ns_2d$X_ns[,2],y=data_ns_2d$Y1_ns)

save(data_ns_2a,file="datans2_1k5.Rdata")
save(data_ns_2b,file="datans2_1k20.Rdata")
save(data_ns_2c,file="datans2_5k5.Rdata")
save(data_ns_2d,file="datans2_5k20.Rdata")


# generating network (spatial) data with interference

##  spatial network 
net5k<-igraph::sample_grg(5000,radius=0.03,torus=FALSE,coords=TRUE)
net1k<-igraph::sample_grg(1000,radius=0.05,torus=FALSE,coords=TRUE)
# check that there's no disconnected nodes.
plot(net1k,vertex.size = 2,vertex.label = NA,
     axes=TRUE,margin=c(0,0,0,0),main="",xlim=c(-1,1),ylim=c(-1,1))
hist(degree(net1k)) 
net5k_adj<-as_adjacency_matrix(net5k)
net1k_adj<-as_adjacency_matrix(net1k)
# save(net5k_adj,file="net5k_adj.Rdata")
# save(net1k_adj,file="net1k_adj.Rdata")
# save(net5k,file="net5k_structure.Rdata")
# save(net1k,file="net1k_structure.Rdata")

load("net5k_structure.Rdata")
load("net1k_structure.Rdata")
load("net5k_adj.Rdata")
load("net1k_adj.Rdata")



##  net1. net2
# spillover of outcome. 
# do not run
gen_nt1_spillY <- function(net,n,d,data_ns_1){
  if(is.null(data_ns_1)){
    data_ns_1<-gen_ns_1(n,d)
  } else {data_ns_1 <- data_ns_1}
  g1 <- igraph::as_data_frame(net,what="vertices")
  colnames(g1)<-c("coordx","coordy")
  g1$Z <- data_ns_1$T_ns_1
  g1$YZ <-data_ns_1$Y_ns_1
  g1$G <- 0
  g1$YG <- 0
  
  for (i in 1:n){
    neigh = igraph::neighbors(net1,v = i)
    g1$G[i] =  mean(g1$Z[neigh]) # can also use inverse dist w
    g1$YG[i] = mean(g1$YZ[neigh]) # spillover of outocme YG
  }
  YG1_nt_1<-g1$YG  
  YG0_nt_1<-0
  YZ1_nt_1<-Y1_ns_1 
  YZ0_nt_1<-Y0_ns_1
  
  Y00_nt_1<-YZ0_nt_1+YG0_nt_1
  Y01_nt_1<-YZ0_nt_1+YG1_nt_1
  Y10_nt_1<-YZ1_nt_1+YG0_nt_1
  Y11_nt_1<-YZ1_nt_1+YG1_nt_1
  
  A_adj<-as_adjacency_matrix(net)
  A<-as.matrix(A_adj)
  Xnei_nt_1 <- get_Xnei(data_ns_1$X_ns_1,A,"adj")
  ps_nt_1g <- pnorm(Xnei_nt_1[,1]) 
  ps_nt_1z <- ps_ns_1
  ps_nt_1sum <- data.frame(z=ps_nt_1z,g=ps_nt_1g)
  G_nt_1 <- sapply(ps_nt_1g,function(x){rbinom(1,1,prob=x)})
  g1$YG <- ifelse(G_nt_1==1, g1$YG, 0)
  # summarise 
  Y_nt_1<-g1$YZ+g1$YG 
  X_nt_1<-data_ns_1$X_ns_1
  Z_nt_1<-g1$Z # Z_nt_1 <- sapply(ps_ns_1,function(x){rbinom(1,1,prob=x)})
  T_nt_1<-paste0("t",Z_nt_1,G_nt_1)
  
  return(list(X_nt_1=X_nt_1,Y_nt_1=Y_nt_1,Y00_nt_1=Y00_nt_1,Y01_nt_1=Y01_nt_1,
              Y10_nt_1=Y10_nt_1,Y11_nt_1=Y11_nt_1,
              T_nt_1=T_nt_1,ps_nt_1sum=ps_nt_1sum))
}



# spillover of treatment
gen_nt1 <- function(net,n,d,data_ns,linear=TRUE){
  
  if(is.null(data_ns)){
    if(linear==TRUE){
      data_ns<-gen_ns(n,d,linear=TRUE)
    } else {data_ns<-gen_ns(n,d,linear=FALSE)}
  } else {data_ns <- data_ns}
  
  YZ1<-data_ns$Y1_ns 
  YZ0<-data_ns$Y0_ns
  YG1<-0.5*YZ1 # outcome from spillover of treatment
  YG0<-0
  
  Y00<-YZ0+YG0
  Y01<-YZ0+YG1
  Y10<-YZ1+YG0
  Y11<-YZ1+YG1
  
  A_adj<-as_adjacency_matrix(net)
  A<-as.matrix(A_adj)
  Xnei_nt <- get_Xnei(data_ns$X_ns,A,"adj")
  ps_nt_1g <- pnorm(Xnei_nt[,1]) 
  ps_nt_1z <- data_ns$ps_ns
  ps_nt <- data.frame(z=ps_nt_1z,g=ps_nt_1g)
  G_nt <- sapply(ps_nt_1g,function(x){rbinom(1,1,prob=x)})
  YG <- ifelse(G_nt==1,YG1, 0)
  # summarise and save data
  Y_nt<-data_ns$Y_ns+YG 
  X_nt<-data_ns$X_ns
  Z_nt<-data_ns$T_ns 
  T_nt<-paste0("t",Z_nt,G_nt)
  
  return(list(X_nt=X_nt,Y_nt=Y_nt,Y00=Y00,Y01=Y01,
              Y10=Y10,Y11=Y11,Z_nt=Z_nt,
              T_nt=T_nt,ps_nt=ps_nt))
}

data_nt_1a <- gen_nt1(net1k,1000,5,data_ns_1a,linear=TRUE)
data_nt_1b <- gen_nt1(net1k,1000,20,data_ns_1b,linear=TRUE)
data_nt_1c <- gen_nt1(net5k,5000,5,data_ns_1c,linear=TRUE)
data_nt_1d <- gen_nt1(net5k,5000,20,data_ns_1d,linear=TRUE)

hist(data_nt_1d$X_nt[,1])
hist(data_nt_1d$ps_nt[,1])
hist(data_nt_1d$Y_nt)
sd(data_nt_1d$Y_nt) # 7.34
plot(x=data_nt_1d$X_nt[,1],y=data_nt_1d$Y11)

save(data_nt_1a,file="datant1_1k5.Rdata")
save(data_nt_1b,file="datant1_1k20.Rdata")
save(data_nt_1c,file="datant1_5k5.Rdata")
save(data_nt_1d,file="datant1_5k20.Rdata")

data_nt_2a <- gen_nt1(net1k,1000,5,data_ns_2a,linear=FALSE)
data_nt_2b <- gen_nt1(net1k,1000,20,data_ns_2b,linear=FALSE)
data_nt_2c <- gen_nt1(net5k,5000,5,data_ns_2c,linear=FALSE)
data_nt_2d <- gen_nt1(net5k,5000,20,data_ns_2d,linear=FALSE)

hist(data_nt_2d$X_nt[,1])
hist(data_nt_2d$ps_nt[,1])
hist(data_nt_2d$Y_nt)
sd(data_nt_2d$Y_nt) # 12.08
plot(x=data_nt_2d$X_nt[,1],y=data_nt_2d$Y11)

save(data_nt_2a,file="datant2_1k5.Rdata")
save(data_nt_2b,file="datant2_1k20.Rdata")
save(data_nt_2c,file="datant2_5k5.Rdata")
save(data_nt_2d,file="datant2_5k20.Rdata")



## net3. 

# spatial trend
set.seed(123)
grid <- expand.grid(1:100, 1:100)
names(grid)<-c("x","y")
g.dummy1 <- gstat(formula=z~1+x+y, locations=~x+y, dummy=T, beta=c(1,0.01,0.005),
                  model=vgm(psill=0.025, range=15, model='Exp'), nmax=20)
coeff_grid <- predict(g.dummy1, newdata=grid, nsim=4)
sp::gridded(coeff_grid) = ~x+y
# save(coeff_grid,file="coeff_grid.Rdata") 


#load("coeff_grid.Rdata")

# net3 data
gen_nt3 <- function(net,n,d,data_ns,coeff_grid){
  g1 <- igraph::as_data_frame(net,what="vertices")
  colnames(g1)<-c("coordx","coordy")
  sim<-data.frame(coeff_grid$sim1,coeff_grid$sim2,coeff_grid@coords)
  colnames(sim)<-c("coeff1","coeff2","x","y")
  sim$xy<-paste0(sim$x,"-",sim$y)
  coeff<-data.frame(floor(g1$coordx*100)+1,floor(g1$coordy*100)+1)
  colnames(coeff)<-c("x","y")
  coeff$xy<-paste0(coeff$x,"-",coeff$y)
  coeff<-left_join(coeff,sim[,c(1,2,5)],by="xy")
  
  beta_nt <- rep(1, d)
  genY_nt_3 <- function(X){
    10/(1 + exp(coeff$coeff1*(X[, 1]-0.5)))*2/(1 + exp(-1 *(X[, 2] - 0.5)))+
      X%*% beta_nt-10+rnorm(n)
  }
  YZ1 <- genY_nt_3(X=data_ns$X_ns) 
  YZ0<-data_ns$Y0_ns
  YG1<-0.5*YZ1
  YG0<-0
  Y00<-YZ0+YG0
  Y01<-YZ0+YG1
  Y10<-YZ1+YG0
  Y11<-YZ1+YG1
  
  X_nt<-data_ns$X_ns
  A_adj<-as_adjacency_matrix(net)
  A<-as.matrix(A_adj)
  Xnei_nt <- get_Xnei(X_nt,A,"adj")
  ps_g <- pnorm(Xnei_nt[,1]) # closely centered around 0.5
  G_nt <- sapply(ps_g,function(x){rbinom(1,1,prob=x)})
  
  YG <- ifelse(G_nt==1, YG1, 0)
  YZ <- ifelse(data_ns$T_ns == 1, YZ1, YZ0)
  Y_nt<- YZ + YG
  
  Z_nt<-data_ns$T_ns
  T_nt<-paste0("t",Z_nt,G_nt)
  ps_z <- data_ns$ps_ns
  ps_nt <- data.frame(z=ps_z,g=ps_g)
  
  return(list(X_nt=X_nt,Y_nt=Y_nt,Y00=Y00,Y01=Y01,
              Y10=Y10,Y11=Y11,Z_nt=Z_nt,
              T_nt=T_nt,ps_nt=ps_nt))
}

data_nt_3a <- gen_nt3(net1k,1000,5,data_ns_2a,coeff_grid)
data_nt_3b <- gen_nt3(net1k,1000,20,data_ns_2b,coeff_grid)
data_nt_3c <- gen_nt3(net5k,5000,5,data_ns_2c,coeff_grid)
data_nt_3d <- gen_nt3(net5k,5000,20,data_ns_2d,coeff_grid)

hist(data_nt_3c$X_nt[,1])
hist(data_nt_3c$ps_nt[,1])
hist(data_nt_3c$Y_nt)
sd(data_nt_3c$Y_nt) # 8.81
plot(x=data_nt_3c$X_nt[,1],y=data_nt_3c$Y11)

save(data_nt_3a,file="datant3_1k5.Rdata")
save(data_nt_3b,file="datant3_1k20.Rdata")
save(data_nt_3c,file="datant3_5k5.Rdata")
save(data_nt_3d,file="datant3_5k20.Rdata")



## net4. With network autocorrelation. 

# net4 data
gen_nt4 <- function(net,n,d,data_ns){
  
  X_nt <- data_ns$X_ns 
  A_adj<-as_adjacency_matrix(net)
  A<-as.matrix(A_adj)
  Xnei_nt <- get_Xnei(X_nt,A,"adj")
  X_nt[,1]<- 0.5*X_nt[,1]+Xnei_nt[,1]
  X_nt[,2]<- 0.5*X_nt[,1]+Xnei_nt[,2] 
  
  beta_nt <- rep(1, d)
  genY_nt_4 <- function(X){
    (10/(1+exp(-12*(X[, 1]-0.5))))*(2/(1+exp(-12*(X[, 2]-0.5))))+
      X%*%beta_nt-10+rnorm(n)
  }
  YZ1 <- genY_nt_4(X_nt)
  YZ0 <- rnorm(n)
  ps_z <- pnorm(X_nt[, 1]+X_nt[, 2])
  Z_nt <- sapply(ps_z,function(x){rbinom(1,1,prob=x)})
  YZ <- ifelse(Z_nt==1, YZ1, YZ0)
  
  YG1<-0.5*YZ1  
  YG0<-0 
  Y00<-YZ0+YG0
  Y01<-YZ0+YG1
  Y10<-YZ1+YG0
  Y11<-YZ1+YG1
  
  Xnei_nt <- get_Xnei(X_nt,A,"adj") # update Xnei
  ps_g <- pnorm(Xnei_nt[,1])
  G_nt <- sapply(ps_g,function(x){rbinom(1,1,prob=x)})
  YG <- ifelse(G_nt==1, YG1, 0)
  Y_nt<-YZ+YG 
  
  T_nt<-paste0("t",Z_nt,G_nt)
  ps_nt <- data.frame(z=ps_z,g=ps_g)
  
  return(list(X_nt=X_nt,Y_nt=Y_nt,Y00=Y00,Y01=Y01,
              Y10=Y10,Y11=Y11,Z_nt=Z_nt,
              T_nt=T_nt,ps_nt=ps_nt))  
}

data_nt_4a <- gen_nt4(net1k,1000,5,data_ns_2a)
data_nt_4b <- gen_nt4(net1k,1000,20,data_ns_2b)
data_nt_4c <- gen_nt4(net5k,5000,5,data_ns_2c)
data_nt_4d <- gen_nt4(net5k,5000,20,data_ns_2d)

hist(data_nt_4d$X_nt[,1])
hist(data_nt_4d$ps_nt[,1])
hist(data_nt_4d$Y_nt)
sd(data_nt_4d$Y_nt) # 9.68
plot(x=data_nt_4d$X_nt[,1],y=data_nt_4d$Y11)

save(data_nt_4a,file="datant4_1k5.Rdata")
save(data_nt_4b,file="datant4_1k20.Rdata")
save(data_nt_4c,file="datant4_5k5.Rdata")
save(data_nt_4d,file="datant4_5k20.Rdata")



## plot data

# plotting net3
# the spatial fields
spplot(coeff_grid)

# X, Y, Z and G
mapcheck<-cbind.data.frame(x=coeff$x,y=coeff$y,X_ns_1,Z_nt_3,G_nt_3,
                           Y_nt_3,ps_nt_3g,ps_nt_2z)

ggplot(mapcheck,aes(x=x,y=y,color=`1`))+
geom_point(size=4)+xlim(0,50)+ylim(0,50)+
scale_colour_gradient(low = "darkgreen", high = "yellow",name="Covariate X1")+
theme_minimal()+
labs(x="coordinate 1",y="coordinate 2",title="Spatially random X")

ggplot(mapcheck,aes(x=x,y=y,color=Y_nt_3))+
geom_point(size=4)+xlim(0,50)+ylim(0,50)+
scale_colour_gradient(low = "darkgreen", high = "yellow",name="Y outcome")+
theme_minimal()+
labs(x="coordinate 1",y="coordinate 2",title="Outcome from patially random X")

ggplot(mapcheck,aes(x=x,y=y,color=Z_nt_3))+
  geom_point(size=4)+xlim(0,50)+ylim(0,50)+
  scale_colour_gradient(low = "darkgreen", high = "yellow",name="Z level")+
  theme_minimal()+
  labs(x="coordinate 1",y="coordinate 2",title="Z treatment from spatially random X")

ggplot(mapcheck,aes(x=x,y=y,color=G_nt_3))+
  geom_point(size=4)+xlim(0,50)+ylim(0,50)+
  scale_colour_gradient(low = "darkgreen", high = "yellow",name="G level")+
  theme_minimal()+
  labs(x="coordinate 1",y="coordinate 2",title="G treatment from spatially random X")

ggplot(coeff,aes(x=x,y=y,color=coeff1))+
  geom_point(size=3)+
  scale_colour_gradient(low = "darkgreen", high = "yellow")+
  theme_minimal()+
  labs(x="coordinate 1",y="coordinate 2",title="Spatially varying causal coefficient")



# plotting net4
mapcheck2<-cbind.data.frame(x=coeff$x,y=coeff$y,Xnei_nt_4,X_nt_4,Z_nt_4,G_nt_4,Y_nt_4,
                            ps_nt_4g,ps_nt_4z)

ggplot(mapcheck2,aes(x=x,y=y,color=1.4*`V1`))+
  geom_point(size=4)+xlim(0,50)+ylim(0,50)+
  scale_colour_gradient(low = "darkgreen", high = "yellow",name="Smoothed X1")+
  theme_minimal()+labs(x="coordinate 1",y="coordinate 2",title="Spatially smoothed X")

ggplot(mapcheck2,aes(x=x,y=y,color=Z_nt_4))+
  geom_point(size=4)+xlim(0,50)+ylim(0,50)+
  scale_colour_gradient(low = "darkgreen", high = "yellow",name="Z level")+
  theme_minimal()+
  labs(x="coordinate 1",y="coordinate 2",title="Z treatment from spatially smoothed X")

ggplot(mapcheck2,aes(x=x,y=y,color=G_nt_4))+
  geom_point(size=4)+xlim(0,50)+ylim(0,50)+
  scale_colour_gradient(low = "darkgreen", high = "yellow",name="G level")+
  theme_minimal()+
  labs(x="coordinate 1",y="coordinate 2",title="G treatment from spatially smoothed X")

ggplot(mapcheck2,aes(x=x,y=y,color=Y_nt_4))+
  geom_point(size=4)+xlim(0,50)+ylim(0,50)+
  scale_colour_gradient(low = "darkgreen", high = "yellow",name="Y outcome")+
  theme_minimal()+
  labs(x="coordinate 1",y="coordinate 2",title="Outcome from spatially smoothed X")

save(ps_nt_4g,ps_nt_4z,Y00_nt_4,Y01_nt_4,Y10_nt_4,Y11_nt_4,Y_nt_4,X_nt_4,Z_nt_4,G_nt_4,T_nt_4,file="datanet4new_5k.Rdata")


