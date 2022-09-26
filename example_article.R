library(mexhaz)
library(spBayesSurv)
library(survival)
data(LeukSurv)
LeukSurv$timeY <- (LeukSurv$time / 365.25)

nodel <- median(LeukSurv$timeY[LeukSurv$cens==1])
modelBS <- mexhaz(formula = Surv(time=timeY , event=cens) ~ age + sex + wbc + tpi + nph(age),data=LeukSurv, base = "exp.bs", degree = 3, knots=nodel )

#Construction of the model predicted matrix
survieBS <- vector()
predBS <- vector()
for (i in 1:nrow(LeukSurv)) {
  predBS <- predict(modelBS , time.pts= sort(LeukSurv$timeY),  
                    data.val = data.frame(age = LeukSurv$age[i],sex=LeukSurv$sex[i], wbc=LeukSurv$wbc[i],tpi=LeukSurv$tpi[i]))$results$surv
  survieBS <-cbind(survieBS ,predBS )
}
TsurvieBS <- sort(LeukSurv$timeY)

#Construction of the censoring matrix
coxcens <- coxph(Surv(timeY, 1-cens) ~ 1,data = LeukSurv) 
pgroupe <- ncol(coxcens[["y"]])
tcens <- ycens[, pgroupe - 1]
dcens <- ycens[, pgroupe]
xcens  <- coxcens$linear.predictors
coxcens  <- coxph(Surv(tcens , dcens ) ~ xcens)             
sfcens  <- survfit(coxcens, newdata = data.frame(x = xcens ), type = "kalbfl")
tcens  <- sfcens$time
censcox <- sfcens$surv

utimes <- sort(unique(LeukSurv$timeY))


BrierScore <- function(times, status, timesurv, matSurv, timecens, matCens, tseq, N){
  calcBS <- vector()
  predcentime <- vector()
  predcent <- vector()
  j<-1
  for(t in tseq){
    Itimeinft<- rep(NA,N) # indicator I(time_i <= t and status_i = 1)
    Itimesupt<- rep(NA,N) # indicator I(time_i > t)
    Itimeinft <- ifelse(times<=t & status==1,1,0)
    Itimesupt <- ifelse(times>t,1,0)
    
    indsurv<-tail(which(timesurv<=t),1) # to get the right row in the matrix
    pred<-as.numeric(matSurv[indsurv,]) # get the predictions of the survival matrix
    
    indcens<-tail(which(timecens<=t),1) # to get the right row in the matrix
    
    for (i in 1:N){
      predcentime[i]<-as.numeric(matCens[tail(which(timecens<=times[i]),1),i])[1]# G(Ti)
      ifelse(predcentime[i]==0,predcentime[i]<- -Inf,predcentime[i])
      ifelse(is.na(predcentime[i]),predcentime[i]<- -Inf,predcentime[i])
      indice <- tail(which(timecens<=timesurv[i]),1)
      predcent[i]<-as.numeric(matCens[tail(which(timecens<=t),1),indice])[1]# G(t)
    }
    #get the times of the group
    ifelse(predcent==0,predcent<- -Inf,predcent)
    timeinft <- Itimeinft*(1/predcentime)*((0-pred)^2)
    timesupt <-Itimesupt*(1/predcent)*((1-pred)^2)
    
    calcBS[j]<-mean(timeinft + timesupt )
    j <- j+1
  }
  return(data.frame(tseq,calcBS))
}



briersBS <-BrierScore(LeukSurv$timeY, LeukSurv$cens, TsurvieBS , survieBS , tcens, censcox, utimes ,N=length(TsurvieBS))
plot(x=briersBS[,1],y=briersBS[,2])
# Graphique vide
plot(0,0,col="white",xlab="Time",ylab="Prediction error",xlim=c(0,8),ylim=c(0,0.3), cex.lab=1.5,cex.axis=1.5,
     cex.main=1.5,main="Prediction error curves")
lines(briersBS[,1],briersBS[,2],col="black",type="s",lty=2)



calibQuantile <- function(quantiles, vecpred, vecobs){
  dfCalibration <- data.frame(row.names = c(1:length(quantiles)))
  for(i in 1:length(quantiles)){
    if(i==1){
      indice <- as.numeric(names(vecpred[which(vecpred <= quantiles[i])]))
      dfCalibration$observed[i] <- mean(vecobs[indice],na.rm=TRUE)
      dfCalibration$predict[i] <- mean(vecpred[which(vecpred <= quantiles[i])])
    }
    else{
      indice <- as.numeric(names(vecpred[which(vecpred > quantiles[i-1] & vecpred <= quantiles[i])]))
      dfCalibration$observed[i] <- mean(vecobs[indice],na.rm=TRUE)
      dfCalibration$predict[i] <- mean(vecpred[which(vecpred > quantiles[i-1] & vecpred <= quantiles[i])])
    }   
  }
  return(dfCalibration)
}


#Calibration at time 3 years
LeukSurv$TH3 <- ifelse(LeukSurv$timeY<=3 & LeukSurv$cens==1,1,
                       ifelse(LeukSurv$timeY<=3 & LeukSurv$cens==0,NA,0))
indiceBS<-which(TsurvieBS==min(TsurvieBS[TsurvieBS>=3]))[1]
pred3BS<-1-survieBS[indiceBS,]  
names(pred3BS)<-c(1:length(pred3BS))
pred3BS_sort <-sort(pred3BS)
quantileBS<-quantile(pred3BS_sort,probs = seq(0.1,1,0.1))
calibBS<-calibQuantile(quantileBS,vecpred=pred3BS_sort,
                       vecobs=LeukSurv$TH3)

plot(0,0,col="white",xlim=c(0,8),ylim=c(0,1),xlab = "Time",
     ylab="Survival",main="Survival curves for different\nvalue of Sex",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
levSex<-sort(unique(LeukSurv$sex)) # 0 : "Female" - 1 : "Male"
for(L in levSex){
  indSex <- which(LeukSurv$sex==L)
  survBSSex<-vector()
  for(i in 1:length(sort(unique(LeukSurv$timeY)))){
    survBSSex[i] <- mean(sfBS[i,indSex])        #BS model
  }
  lines(utimes ,survBSSex,col=L+4,lty=3,lwd=2)
  
  
  # Kaplan Meier
  KMSex <- LeukSurv[LeukSurv$sex==L,]
  lines(survfit(Surv(KMSex$timeY,KMSex$cens)~1),data=KMSex,conf.int=FALSE,col="gray40",lty=1,lwd=2)
}



AUCt.TPFP <-function (marker, Stime, status, predict.time) 
{
  eta <- marker
  Target <- predict.time
  at.risk <- ((Stime >= Target))
  the.eta <- eta[at.risk] #we keep only patients still at risk at time t 
  n <- length(the.eta)  
  the.dead <- (Stime == Target) & (status == 1) # patient who are deceaded at time t 
  the.dead <- the.dead[at.risk] # 
  n.dead <- sum(the.dead)
  
  
  ooo <- order(the.eta)
  eta.out <- the.eta[ooo]
  the.dead.out <- the.dead[ooo] # vector of deceaded patients among those who still at risk at time t 
  the.dead.out_inv <- 1 - the.dead.out
  TP = NULL
  FP = NULL
  # c<-eta.out[the.dead.out]
  for (c in eta.out){
    TP = c(TP,sum((eta.out>c)*exp(eta.out))/sum(exp(eta.out)))
    FP = c(FP,sum((eta.out>c)/(n-n.dead)* the.dead.out_inv)) # 
  }
  length(FP)
  TP <- c(1, TP, 0)
  FP <- c(1, FP, 0)
  # Area under the curve of false positive and true positive ;
  dFP <- abs(FP[-1] - FP[-length(FP)])
  aTP <- 0.5 * (TP[-1] + TP[-length(TP)])
  area <- sum(dFP * aTP)
  out <- list(marker = eta.out, TP = TP, FP = FP, AUC = area)
  return(out)
} 



AUCBS = NULL 
for (i in 1: NROW(utimes)) {
  bashaz = predict(modelBS, time.pts = utimes[i], 
                   data.val = data.frame("age"=0 ,"sex"=0 , "wbc"=0,  "tpi"=0), data=LeukSurv)$results$hazard
  haz1 = predict(modelBS, time.pts =  utimes[i], data.val = 
                   data.frame(LeukSurv[, c("age", "sex", "wbc","tpi")])) $results$hazard
  new.eta1 = log(haz1/bashaz)
  
  out1 = MyCoxWeights(marker = new.eta1, Stime = LeukSurv$timeY, 
                      status = LeukSurv$cens,predict.time = utimes[i])
  
  AUCBS = c(AUCBS, out1$AUC)
}

plot(utimes,AUCBS,type="l",col="gray",lwd=2,lty=1,ylim = c(0.6,0.8),xlim = c(0,10)) 

  
