#Load Packages
##############################

library(MASS)
library(psych)
library(EnvStats)
library(data.table)
library(psych)
library(glmmLasso)
library(lme4)

#Function AR(1)
##############################

c.mat<-function(phi,ord){
  c=matrix(rep(phi, ord*ord), nrow = ord, ncol = ord)
  power <- abs(outer(1:ord, 1:ord, "-"))
  c<- c^power
  return(c)
}

#Function for MtLMM 
##############################

MtLMM = function(id, n, si, r, ni, tik, ti, xi1, xi2, xi, zi1, zi2, zi, yi1, yi2, yi)
{

#E-step of ECM 
##############################

 phihat0 <- 0
  
  Omegahati0<-list()
  for(i in 1:n){Omegahati0[[i]]=c.mat(phihat0,si[[i]])}

  Sigmahat0 <- diag(1,r)
  
  Rihat0 <- list()
  for(i in 1:n){Rihat0[[i]] = Sigmahat0%x%Omegahati0[[i]]}
  
  Dhat0 <- diag(1,dim(zi[[1]])[2])
  
  Lambdahati0 <- list()
  for(i in 1:n){Lambdahati0[[i]] = zi[[i]]%*%Dhat0%*%t(zi[[i]])+Rihat0[[i]]}
  
  xxi <- list()
  for(i in 1:n){xxi[[i]] = t(xi[[i]])%*%xi[[i]]}
  
  xyi <- list()
  for(i in 1:n){xyi[[i]] = t(xi[[i]])%*% yi[[i]]}

  betahat0 <- ginv(Reduce('+', xxi))%*%(Reduce('+', xyi))
  
  Deltahati0 <- list()
  for(i in 1:n){Deltahati0[[i]] = t(yi[[i]]-xi[[i]]%*%betahat0)%*%ginv(Lambdahati0[[i]])%*%(yi[[i]]-xi[[i]]%*%betahat0)}
  
  nuhat0 <- 30
  
  tauhati0 <- list()
  for(i in 1:n){tauhati0[[i]]=(nuhat0+ni[[i]])/(nuhat0+Deltahati0[[i]])}

  kapahati0 <- list()
  for(i in 1:n){kapahati0[[i]] <- digamma((nuhat0+ni[[i]])/2)-log((nuhat0+Deltahati0[[i]])/2)}
  
  bhati0 <- list()
  for(i in 1:n){bhati0[[i]] = Dhat0%*%t(zi[[i]])%*%ginv(Lambdahati0[[i]])%*%(yi[[i]]-xi[[i]]%*%betahat0)}
  
  bhati01 <- list()
  for(i in 1:n){bhati01[[i]] = bhati0[[i]][1:2,]}
  
  bhati02 <- list()
  for(i in 1:n){bhati02[[i]] = bhati0[[i]][3:4,]}
  
  Vhatbi0 <- list()
  for(i in 1:n){Vhatbi0[[i]] = ginv(ginv(Dhat0)+t(zi[[i]])%*%ginv(Rihat0[[i]])%*%zi[[i]])}
  
  Vhatbi011 <- list()
  for(i in 1:n){Vhatbi011[[i]] = Vhatbi0[[i]][1:2,1:2]}
  
  Vhatbi012 <- list()
  for(i in 1:n){Vhatbi012[[i]] = Vhatbi0[[i]][1:2,3:4]}
  
  Vhatbi021 <- list()
  for(i in 1:n){Vhatbi021[[i]] = Vhatbi0[[i]][3:4,1:2]}
  
  Vhatbi022 <- list()
  for(i in 1:n){Vhatbi022[[i]] = Vhatbi0[[i]][3:4,3:4]}
  
  Bhati0 <- list()
  for(i in 1:n){Bhati0[[i]] = bhati0[[i]]%*%t(bhati0[[i]])+Vhatbi0[[i]]}  
  
  Psihati0 <- list()
  for(i in 1:n){Psihati0[[i]]=(yi[[i]]-xi[[i]]%*%betahat0-zi[[i]]%*%bhati0[[i]])%*%t(yi[[i]]-xi[[i]]%*%betahat0-zi[[i]]%*%bhati0[[i]])+zi[[i]]%*%Vhatbi0[[i]]%*%t(zi[[i]])}
  
  betahattek <- matrix(NA,ncol=length(betahat0),nrow=tek)
  betahat <- c(rep(NA,length(betahat0)))
  phihat <- c(rep(NA,length(phihat0)))
  Sigmahat <- matrix(NA,ncol=r,nrow=r)
  Dhat <- matrix(NA,ncol=r,nrow=r)
  nuhat <- c(rep(NA,length(nuhat0)))

#LOOP
###########################
 
 j <- 1
  while (j <= tek) {
    
    tauxRx <- list()
    for(i in 1:n){tauxRx[[i]]<-as.numeric(tauhati0[[i]])*t(xi[[i]])%*%ginv(Rihat0[[i]])%*%xi[[i]]}
    
    tauxRy_z <-list()
    for(i in 1:n){tauxRy_z[[i]]<-as.numeric(tauhati0[[i]])*t(xi[[i]])%*%ginv(Rihat0[[i]])%*%(yi[[i]]-zi[[i]]%*%bhati0[[i]])}
    
    betahat <- ginv(as.matrix(Reduce('+', tauxRx)))%*%(Reduce('+', tauxRy_z))
    
    betahat1 <- betahat[1:dim(xi1[[1]])[2]]
    betahat2 <- betahat[(dim(xi1[[1]])[2]+1):dim(xi[[1]])[2]]
    
    betahattek[j,] <- betahat
    
    Dhat <- Reduce('+', Bhati0)/n
    
    eijhat <- list()
    for(i in 1:n){eijhat[[i]] = yi[[i]]-xi[[i]]%*%betahat-zi[[i]]%*%bhati0[[i]]}
    
    eijhat1 <- list()
    for(i in 1:n){eijhat1[[i]] = yi1[[i]]-xi1[[i]]%*%betahat1-zi1[[i]]%*%bhati01[[i]]}
    
    eijhat2 <- list()
    for(i in 1:n){eijhat2[[i]] = yi2[[i]]-xi2[[i]]%*%betahat2-zi2[[i]]%*%bhati02[[i]]}
    
    psihati11 <- list()
    for(i in 1:n){psihati11[[i]] = eijhat1[[i]]%*%t(eijhat1[[i]])-zi1[[i]]%*%Vhatbi011[[i]]%*%t(zi1[[i]])}
    
    psihati12 <- list()
    for(i in 1:n){psihati12[[i]] = eijhat1[[i]]%*%t(eijhat2[[i]])-zi1[[i]]%*%Vhatbi012[[i]]%*%t(zi2[[i]])}
    
    psihati21 <- list()
    for(i in 1:n){psihati21[[i]] = eijhat2[[i]]%*%t(eijhat1[[i]])-zi2[[i]]%*%Vhatbi021[[i]]%*%t(zi1[[i]])}
    
    psihati22 <- list()
    for(i in 1:n){psihati22[[i]] = eijhat2[[i]]%*%t(eijhat2[[i]])-zi2[[i]]%*%Vhatbi022[[i]]%*%t(zi2[[i]])}
    
    ompshati11 <- list()
    for(i in 1:n){ompshati11[[i]] = tr(ginv(Omegahati0[[i]])%*%psihati11[[i]])}
    
    ompshati22=list()
    for(i in 1:n){ompshati22[[i]]=tr(ginv(Omegahati0[[i]])%*%psihati22[[i]])}
    
    ompshati12 <- list()
    for(i in 1:n){ompshati12[[i]] = tr(ginv(Omegahati0[[i]])%*%(psihati12[[i]]+psihati21[[i]]))}
    
    sigma11 <- (1/Reduce('+', si))*Reduce('+', ompshati11)
    
    sigma22 <- (1/Reduce('+', si))*Reduce('+', ompshati22)
    
    sigma12 <- (1/(2*Reduce('+', si)))*Reduce('+', ompshati12)
    
    Sigmahat <- matrix(c(sigma11, sigma12, sigma12, sigma22),ncol=2,nrow=2)
    
    Psihati.2 <- list()
    for(i in 1:n){Psihati.2[[i]] = (yi[[i]]-xi[[i]]%*%betahat-zi[[i]]%*%bhati0[[i]])%*%t(yi[[i]]-xi[[i]]%*%betahat-zi[[i]]%*%bhati0[[i]])+zi[[i]]%*%Vhatbi0[[i]]%*%t(zi[[i]])}
    
    phihat <- 1/exp(1)
    
    f2<-function(nu){
      li<-vector()
      for(i in 1:n){li[i]=-nu*log(nu/2)+2*lgamma(nu/2)-nu*(kapahati0[[i]]-tauhati0[[i]])}
      L=sum(li)
      return(L)
    }
    
    nuhat=nlminb(start=4,objective=f2,lower=1, upper=Inf)$par

    Omegahati<-list()
    for(i in 1:n){Omegahati[[i]]=c.mat(phihat,si[[i]])}

    Rihat <- list()
    for(i in 1:n){Rihat[[i]] = Sigmahat%x%Omegahati[[i]]}
    
    Lambdahati <- list() 
    for(i in 1:n){Lambdahati[[i]] = zi[[i]]%*%Dhat%*%t(zi[[i]])+Rihat[[i]]}
    
    Deltahati <- list()
    for(i in 1:n){Deltahati[[i]] = t(yi[[i]]-xi[[i]]%*%betahat)%*%ginv(Lambdahati[[i]])%*%(yi[[i]]-xi[[i]]%*%betahat)}
    
    tauhati<-list()
    for(i in 1:n){tauhati[[i]]=(nuhat0+ni[[i]])/(nuhat0+Deltahati0[[i]])}
    
    kapahati=list()
    for(i in 1:n){kapahati[[i]]=digamma((nuhat+ni[[i]])/2)-log((nuhat+Deltahati[[i]])/2)}
        
    bhati <- list()
    for(i in 1:n){bhati[[i]] = Dhat%*%t(zi[[i]])%*%ginv(Lambdahati[[i]])%*%(yi[[i]]-xi[[i]]%*%betahat)}
    
    bhati1 <- list()
    for(i in 1:n){bhati1[[i]] = bhati[[i]][1:2,]}
    
    bhati2 <- list()
    for(i in 1:n){bhati2[[i]] = bhati[[i]][3:4,]}
    
    Vhatbi <- list()
    for(i in 1:n){Vhatbi[[i]] = ginv(ginv(Dhat)+t(zi[[i]])%*%ginv(Rihat[[i]])%*%zi[[i]])}
    
    Vhatbi11 <- list()
    for(i in 1:n){Vhatbi11[[i]] = Vhatbi[[i]][1:2,1:2]}
    
    Vhatbi12 <- list()
    for(i in 1:n){Vhatbi12[[i]] = Vhatbi[[i]][1:2,3:4]}
    
    Vhatbi21 <- list()
    for(i in 1:n){Vhatbi21[[i]] = Vhatbi[[i]][3:4,1:2]}
    
    Vhatbi22 <- list()
    for(i in 1:n){Vhatbi22[[i]] = Vhatbi[[i]][3:4,3:4]}
    
    Bhati <- list()
    for(i in 1:n){Bhati[[i]] = bhati[[i]]%*%t(bhati[[i]])+Vhatbi[[i]]}  
    
    Psihati <- list()
    for(i in 1:n){Psihati[[i]] = (yi[[i]]-xi[[i]]%*%betahat-zi[[i]]%*%bhati[[i]])%*%t(yi[[i]]-xi[[i]]%*%betahat-zi[[i]]%*%bhati[[i]])+zi[[i]]%*%Vhatbi[[i]]%*%t(zi[[i]])}
    #####################
    j <- j + 1
    if (all(abs(betahat-betahat0) < tol)) break
    
    phihat0 <- phihat
    
    for(i in 1:n){Omegahati0[[i]] = Omegahati[[i]]}
    Sigmahat0 <- Sigmahat
    for(i in 1:n){Rihat0[[i]] = Rihat[[i]]}
    Dhat0 <- Dhat
    for(i in 1:n){Lambdahati0[[i]] = Lambdahati[[i]]}
    betahat0 <- betahat
    for(i in 1:n){Deltahati0[[i]] = Deltahati[[i]]}
    nuhat0=nuhat
    for(i in 1:n){tauhati0[[i]]=tauhati[[i]]}
    for(i in 1:n){kapahati0[[i]]=kapahati[[i]]}
    for(i in 1:n){bhati0[[i]] = bhati[[i]]}
    for(i in 1:n){bhati01[[i]] = bhati1[[i]]}
    for(i in 1:n){bhati02[[i]] = bhati2[[i]]}
    for(i in 1:n){Vhatbi0[[i]] = Vhatbi[[i]]}
    for(i in 1:n){Vhatbi011[[i]] = Vhatbi11[[i]]}
    for(i in 1:n){Vhatbi012[[i]] = Vhatbi12[[i]]}
    for(i in 1:n){Vhatbi021[[i]] = Vhatbi21[[i]]}
    for(i in 1:n){Vhatbi022[[i]] = Vhatbi22[[i]]}
    for(i in 1:n){Bhati0[[i]] = Bhati[[i]]}
    for(i in 1:n){Psihati0[[i]] = Psihati[[i]]}
    ##########################
    
  }
  #The end of LOOP
  ##############################
  fit <- list()
  #eatimates of parameters
  fit$beta  <- betahat
  fit$D <- Dhat
  fit$Sigma <- Sigmahat
  fit$nu=nuhat
  fit$b <- bhati
  fit$parhat=c(betahat
               ,Dhat[1,1],Dhat[2,1],Dhat[2,2],Dhat[3,1],Dhat[3,2]
               ,Dhat[3,3],Dhat[4,1],Dhat[4,2],Dhat[4,3],Dhat[4,4]
               ,Sigmahat[1,1],Sigmahat[1,2],Sigmahat[2,2]
               ,nuhat)
  
  #estimates of likelohood
  ##############################
  fit$li <- vector()
  for(i in 1:n){fit$li[i]=0.5*(log(det(ginv(Rihat[[i]])))+log(det(ginv(Dhat)))-tr(ginv(Dhat)%*%Bhati[[i]])-tr(ginv(Rihat[[i]])%*%Psihati[[i]])+nuhat*(log(nuhat/2)+kapahati[[i]]-tauhati[[i]]))-lgamma(nuhat/2)}
  fit$L <- sum(na.omit(fit$li))

  #estimates of AIC and BIC
  ##############################
  fit$mm <- length(fit$parhat)
  fit$AIC <- 2*fit$mm-2*fit$L
  fit$BIC <- fit$mm*log(n)-2*fit$L
  
  fit
}
#end of function MLMM
##############################

#Outpout for AIDS Data
##############################

Data <- fread("C:/Users/HP PAVILION/OneDrive/Desktop/Data4.csv")
tol <- 1E-2
tek <- 10
id = Data$Patid 
id <- as.numeric(factor(Data$Patid,levels=unique(Data$Patid)))
n <- length(unique(id))
si <- list()
for(i in 1:length(unique(id))){si[[i]]=length(which(id==unique(id)[i]))}
r <- 2
ni <- list()
for(i in 1:n){ni[[i]]=r*si[[i]]}
tik <- list()
for(i in 1:n){tik[[i]]=Data$Day[id==i]/7}
ti <- list()
for(i in 1:n){ti[[i]]=cbind(tik[[i]],tik[[i]])}
xi1 = xi2 = list()
for(i in 1:n){xi1[[i]]= xi2[[i]]= matrix(c(rep(1,si[[i]]),Data$Day[id==i]/7,(Data$Day[id==i]/7)^2,Data$rna[id==i]),ncol=4)}
xi <- list()
for(i in 1:n){xi[[i]] = as.matrix(bdiag(xi1[[i]],xi2[[i]]))}
zi1 <- list()
for(i in 1:n){zi1[[i]] = matrix(c(rep(1,si[[i]]),tik[[i]]),ncol=2)}
zi2 <- list()
for(i in 1:n){zi2[[i]] = matrix(c(rep(1,si[[i]]),tik[[i]]),ncol=2)}
zi <- list()
for(i in 1:n){zi[[i]] = as.matrix(bdiag(zi1[[i]],zi2[[i]]))}
yi1 <- list()
for(i in 1:n){yi1[[i]] = Data$log10RNA[id==i]}
yi2 <- list()
for(i in 1:n){yi2[[i]] = Data$CD4.0.5[id==i]}
yi <- list()
for(i in 1:n){yi[[i]] = c(yi1[[i]],yi2[[i]])}

Model.MtLMM.Data2 <- MtLMM (id=id, n=n, si=si, r=r, ni=ni, tik=tik, ti=ti, xi1=xi1, xi2=xi2, xi=xi, zi1=zi1, zi2=zi2, zi=zi, yi1=yi1, yi2=yi2, yi=yi)


# SD and MSPE calculation
########################################

id_index <-  as.numeric(factor(Data$Patid ,levels=unique(Data$Patid )))
n_index <- length(unique(id_index))
si_index <- list()
for(i in 1:n_index){si_index[[i]] = length(which(id_index==unique(id_index)[i]))}

r_index <- 2

ni_index <- list()
for(i in 1:n_index){ni_index[[i]] = r_index*si_index[[i]]}

tik_index <- list()
for(i in 1:n_index){tik_index[[i]] = Data$Day[id_index==i]/7}

ti_index <- list()
for(i in 1:n_index){ti_index[[i]] = cbind(tik_index[[i]],tik_index[[i]])}

xi1_index = xi2_index = list()
for(i in 1:n_index){xi1_index[[i]] = xi2_index[[i]] = matrix(c(rep(1,si_index[[i]]),Data$Day[id_index==i]/7,(Data$Day[id_index==i]/7)^2,Data$rna[id_index==i]),ncol=4)}

xi_index <- list()
for(i in 1:n_index){xi_index[[i]] = as.matrix(bdiag(xi1_index[[i]],xi2_index[[i]]))}

zi1_index <- list()
for(i in 1:n_index){zi1_index[[i]] = matrix(c(rep(1,si_index[[i]]),tik_index[[i]]),ncol=2)}

zi2_index <- list()
for(i in 1:n_index){zi2_index[[i]] = matrix(c(rep(1,si_index[[i]]),tik_index[[i]]),ncol=2)}

zi_index <- list()
for(i in 1:n_index){zi_index[[i]]=as.matrix(bdiag(zi1_index[[i]],zi2_index[[i]]))}

yi1_index <- list()
for(i in 1:n_index){yi1_index[[i]] = Data$log10RNA[id_index==i]}

yi2_index <- list()
for(i in 1:n_index){yi2_index[[i]] = Data$CD4.0.5[id_index==i]}

yi_index <- list()
for(i in 1:n_index){yi_index[[i]] = c(yi1_index[[i]],yi2_index[[i]])}

yhatsimtek <- list()
for(kk in 1:n_index){yhatsimtek[[kk]] = matrix(NA,ncol=si_index[[i]],nrow=1)}

betahatsimtek <- matrix(NA,ncol=8,nrow=n_index)

Sigmasimtek <- list()
for(kk in 1:n_index){Sigmasimtek[[kk]] = matrix(NA,ncol=r_index,nrow=r_index)}

Dsimtek <- list()
for(kk in 1:n_index){Dsimtek[[kk]] = matrix(NA,ncol=r_index,nrow=r_index)}

kk <- 1
while (kk <= n_index) {

  data <- Data[-c(which(Data$Patid==kk)),]

  id <- as.numeric(factor(data$Patid,levels=unique(data$Patid)))
  
  n <- length(unique(id))
  
  si <- list()
  for(i in 1:n){si[[i]] = length(which(id==unique(id)[i]))}
  
  r <- 2
  
  ni <- list()
  for(i in 1:n){ni[[i]] = r*si[[i]]}
  
  tik <- list()
  for(i in 1:n){tik[[i]] = data$Day[id==i]/12}
  
  ti <- list()
  for(i in 1:n){ti[[i]] = cbind(tik[[i]],tik[[i]])}
  
  xi1 = xi2 = list()
  for(i in 1:n){xi1[[i]] = xi2[[i]] = matrix(c(rep(1,si[[i]]),data$Day[id==i]/7,(data$Day[id==i]/7)^2,data$rna[id==i]),ncol=4)}

 
  xi <- list()
  for(i in 1:n){xi[[i]] = as.matrix(bdiag(xi1[[i]],xi2[[i]]))}
  
  zi1 <- list()
  for(i in 1:n){zi1[[i]] = matrix(c(rep(1,si[[i]]),tik[[i]]),ncol=2)}
  
  zi2 <- list()
  for(i in 1:n){zi2[[i]] = matrix(c(rep(1,si[[i]]),tik[[i]]),ncol=2)}
  
  zi <- list()
  for(i in 1:n){zi[[i]] = as.matrix(bdiag(zi1[[i]],zi2[[i]]))}
  
  yi1 <- list()
  for(i in 1:n){yi1[[i]] = data$log10RNA[id==i]}
  
  yi2 <- list()
  for(i in 1:n){yi2[[i]] = data$CD4.0.5[id==i]}
  
  yi <- list()
  for(i in 1:n){yi[[i]] = c(yi1[[i]],yi2[[i]])}
  ##################################################
  cross.Model.MtLMM.Data2= MtLMM (id=id, n=n, si=si, r=r, ni=ni, tik=tik, ti=ti, xi1=xi1, xi2=xi2, xi=xi, zi1=zi1, zi2=zi2, zi=zi, yi1=yi1, yi2=yi2, yi=yi)
  ##################################################
  betahatsimtek[kk,] <- cross.Model.MtLMM.Data2$beta
  Sigmasimtek[[kk]] <- cross.Model.MtLMM.Data2$Sigma
  Dsimtek[[kk]] <- cross.Model.MtLMM.Data2$D
  yhatsimtek[[kk]] <- xi_index[[kk]]%*%cross.Model.MtLMM.Data2$beta+t(t(zi_index[[kk]]))%*%(Reduce("+",cross.Model.MtLMM.Data2$b)/n)
  
  kk <- kk+1
}

mspe <- matrix(NA,ncol=n_index,nrow=1)
for(i in 1:n_index){mspe[i] = t(yi_index[[i]]-yhatsimtek[[i]])%*%(yi_index[[i]]-yhatsimtek[[i]])}

MSPE <- median(as.vector(mspe)+3)

parhat <- list()

for(k in 1:n_index){parhat[[k]] <- cbind(t(betahatsimtek[k,])
                                        ,Dsimtek[[k]][1,1],Dsimtek[[k]][2,1],Dsimtek[[k]][2,2],Dsimtek[[k]][3,1],Dsimtek[[k]][3,2]
                                        ,Dsimtek[[k]][3,3],Dsimtek[[k]][4,1],Dsimtek[[k]][4,2],Dsimtek[[k]][4,3],Dsimtek[[k]][4,4]
                                        ,Sigmasimtek[[k]][1,1],Sigmasimtek[[k]][1,2],Sigmasimtek[[k]][2,2])}

P <- length(parhat[[k]])
mat.par <- matrix(unlist(parhat),ncol=P,byrow = TRUE)

sd <- matrix(NA,ncol=P,nrow=1)
for(i in 1:P){sd[i] = sd(mat.par[,i])}

Model.MtLMM.Data2$AIC
Model.MtLMM.Data2$BIC
Model.MtLMM.Data2$Sigma
MSPE
sd


