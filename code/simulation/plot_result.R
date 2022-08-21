library(tidyverse)
library(ggplot2)
library(viridis)


## load data
load("test6_cov_4.Rdata")
load("test7_cov_4.Rdata")
load("test8_cov_4.Rdata")
load("test9_cov_4.Rdata")
load("test10_cov_4.Rdata")


result<-rbind.data.frame(test6_cov,test7_cov,test8_cov,test9_cov,test10_cov)
result$mean_len2<-rowMeans(result[,5:8])
test<-as.factor(c(rep(6,16),rep(7,16),rep(8,16),rep(9,16),rep(10,16)))
result$test<-test

res1k5<-result[which(result$testID%in%c(61,71,81,91,101)),]
res1k20<-result[which(result$testID%in%c(62,72,82,92,102)),]
res5k5<-result[which(result$testID%in%c(63,73,83,93,103)),]
res5k20<-result[which(result$testID%in%c(64,74,84,94,104)),]


## plot coverage
plot1<-function(data){
  ggplot(data=data,aes(x=exp_cov,y=mean_cov,group=test,color=test))+
    geom_line(size=1.5)+
    geom_point(size=3.5)+
    xlim(0.25,1)+ylim(0.25,1)+
    scale_color_viridis(discrete=TRUE)+
    labs(x="target coverage",y="achieved coverage")+
    annotate("segment",x =0.25,xend =1,y=0.25,yend=1,colour="black",size=0.3,alpha=0.8)+
    theme_minimal()
}

plot1(res1k5)
plot1(res1k20)
plot1(res5k5)
plot1(res5k20)

## plot interval lengths
plot2<-function(data){
  ggplot(data=data,aes(x=exp_meanlen,y=mean_len2,group=test,color=test))+
    geom_line(size=1.5)+
    geom_point(size=3.5)+
    xlim(0,40)+ylim(0,65)+
    annotate("segment",x =0,xend =40,y=0,yend=41,colour="black",size=0.3,alpha=0.8)+
    scale_color_viridis(discrete=TRUE)+
    labs(x="oracle length",y="estimate length")+
    theme_minimal()
}

plot2(res1k5)
plot2(res1k20)
plot2(res5k5)
plot2(res5k20)

## plot coverage disparty by counterfactual treatment levels

# with 5k-20 as illustrative example 
res5k20_l<-res5k20[,-c(5:8)]
res5k20_l<-pivot_longer(res5k20_l,cols=c(1:4),names_to = "Tl",values_to = "cov")

temp<-res5k20[,-c(1:4)]
temp<-pivot_longer(temp,cols=c(1:4),names_to = "Tl",values_to = "len")
res5k20_l$len<-temp$len

res5k20_l7<-res5k20_l[which(res5k20_l$test=="7"),]
res5k20_l8<-res5k20_l[which(res5k20_l$test=="8"),]
res5k20_l9<-res5k20_l[which(res5k20_l$test=="9"),]
res5k20_l10<-res5k20_l[which(res5k20_l$test=="10"),]

plot3<-function(data){
  ggplot(data=data,aes(x=exp_cov,cov,group=Tl,color=Tl))+
    geom_line(size=1.5)+
    geom_point(size=3.5)+
    xlim(0.18,1)+ylim(0.18,1)+
    annotate("segment",x =0.2,xend =1,y=0.2,yend=1,colour="black",size=0.3,alpha=0.8)+
    scale_color_viridis(discrete=TRUE)+
    labs(x="target coverage",y="achieved coverage")+
    theme_minimal()
}

plot3(res5k20_l7)
plot3(res5k20_l8)
plot3(res5k20_l9)
plot3(res5k20_l10)

# with test8 as illustrative example 

res8<-result[which(test=="8"),]
res8_l<-res8[,-c(5:8)]
res8_l<-pivot_longer(res8_l,cols=c(1:4),names_to = "Tl",values_to = "cov")

temp<-res8[,-c(1:4)]
temp<-pivot_longer(temp,cols=c(1:4),names_to = "Tl",values_to = "len")
res8_l$len<-temp$len

res8_la<-res8_l[which(res8_l$testID==81),]
res8_lb<-res8_l[which(res8_l$testID==82),]
res8_lc<-res8_l[which(res8_l$testID==83),]
res8_ld<-res8_l[which(res8_l$testID==84),]

plot3(res8_la)
plot3(res8_lb)
plot3(res8_lc)
plot3(res8_ld)
