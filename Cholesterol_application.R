##required packages
library(skewlmm)
library(ghyp)
library(gsl)

source("em.R")  
source("utils.R")  
source("MI.R")
source("utils_envelopes.R")

################################################################################
## real data analysis: Cholesterol
#################################################################################

##read data
datas<-read.table("chol.txt", header = TRUE, sep = "")
datas$year <- (datas$year-5)/10
datas$cholst<-datas$cholst/100

##prepare data
nj<-matrix(0,200,1)
for (j in 1:200) {
nj[j]=sum(datas$newid==j)
}
y<-datas$cholst
x<-cbind(1,datas$sex,datas$age)
z<- cbind(1,datas$year)
m=dim(nj)[1]
p<-dim(x)[2]
q1<-dim(z)[2]

## Initial values
mod2<- smsn.lmm(datas,formFixed=cholst~sex+age,formRandom = ~1+ year,groupVar="newid",   distr = "st",control = lmmControl(max.iter = 30),   depStruct="CI",pAR=1)
          
betaini<-mod2$est$beta
sigmaeini<- mod2$est$sigma2
D1<-mod2$est$ds
Dini1<- matrix(c(D1[1],D1[2],D1[2],D1[3]), 2, 2)
Dini<-Dini1%*%Dini1
Gamaini<- mod2$est$lambda
eta<- mod2$est$nu  
lambda<- -0.5  # case NIG - initial values to be updated later
phis<-c(lambda,eta)

vj <- matrix(0,nrow=m,ncol=1)  # power 1
ivj <- matrix(0,nrow=m,ncol=1)  # power -1
lvj <- matrix(0,nrow=m,ncol=1)   # log vi

for (j in 1:200 ){
    	y1<-y[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ]
	x1<-matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
	z1<-matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
	med<-x1%*%betaini
        GamaF<-z1%*%Gamaini
	Psi<-(z1)%*%(Dini)%*%t(z1)+diag(nj[j])*sigmaeini
        Psi<-(Psi+t(Psi))/2
        iPsi<-solve(Psi)
	dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
        LambdaA<-solve(solve(Dini)+t(z1)%*%z1/sigmaeini)
	Ai<-(diag(q1)-LambdaA%*%t(z1)%*%z1/sigmaeini)%*%Gamaini
        Bi<-LambdaA%*%t(z1)%*%(y1-med)/sigmaeini
        DeltaF<-eta+t(GamaF)%*%iPsi%*%GamaF
      	vj[j] <- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("x"), check.pars = TRUE)
	ivj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("1/x"), check.pars = TRUE)
    	lvj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("logx"), check.pars = TRUE)  
}

phis <- optim(phis, method = "L-BFGS-B", FCi, lower =c(-10,0.01), upper =c(10,10), vj=vj,ivj=ivj,lvj=lvj,hessian=TRUE)$par

lambdaini<-phis[1] # first update on lambda
etaini<-phis[2] # second update on eta

#get initial values
saida=EM.MGHini(nj,y,x,z,betaini,sigmaeini,Dini, Gamaini, lambdaini, etaini,  type="MGH", MaxIter=100)

betaini1=saida$beta1
sigmaeini1=saida$sigmae
lambdaini1=saida$lambda
etaini1=saida$eta
Dini1=saida$D1
Gamaini1=saida$Gama

############################################
## Fit models
############################################

## Generalized Hyperbolic
EstMGH <- EM.MGH(nj,y,x,z,betaini1,sigmaeini1,Dini1, Gamaini1, lambdaini1, etaini1,  type="MGH", precision=1e-6, MaxIter=1000)

## NIG
EstNIG <- EM.MGH(nj,y,x,z,betaini1,sigmaeini1,Dini1, Gamaini1, lambda=-0.5, etaini1,  type="NIG", precision=1e-6, MaxIter=1000)

# Symmetric
EstMGHsym <- EM.MGH(nj,y,x,z,betaini1,sigmaeini1,Dini1, c(0,0), lambdaini1, etaini1,  type="MGHsym", precision=1e-6, MaxIter=1000)

# SN
EstSn <- smsn.lmm(datas,formFixed=cholst~sex+age,formRandom = ~1+ year,groupVar="newid",   distr = "sn",control = lmmControl(tol=1e-5,max.iter = 1000),   depStruct="CI",pAR=1)     

# ST
EstST <- smsn.lmm(datas,formFixed=cholst~sex+age,formRandom = ~1+ year,groupVar="newid",   distr = "st",control = lmmControl(tol=1e-5,max.iter = 1000),   depStruct="CI",pAR=1)    

#################
## Table 3 
#################
t3 <- matrix(NA, ncol=6,nrow=11)
beta_nam <- paste0("beta",0:2)
D_nam <- paste0("D_1",1:ncol(Dini1))
for(i in 2:nrow(Dini1)) D_nam <- c(D_nam,paste0("D_",i,i:ncol(Dini1)))
names <- c(beta_nam,"sigma_e",D_nam,"lambda","eta","gamma_1","gamma_2")
rownames(t3) <- names
colnames(t3) <- c("GH - Est", "GH - se", "GH sym - Est", "GH sym- se", "NIG - Est", "NIG - se")

#MGH
t3[,1] <- round(c(c(EstMGH$beta1),c(EstMGH$sigmae),EstMGH$D1[upper.tri(EstMGH$D1,diag=TRUE)],EstMGH$lambda,EstMGH$eta,c(EstMGH$Gama)),4)
t3[,2] <- round(EstMGH$ep[c(1:7,10,11,8,9)],4)

#MGH sym
t3[,3] <- round(c(c(EstMGHsym$beta1),c(EstMGHsym$sigmae),EstMGHsym$D1[upper.tri(EstMGHsym$D1,diag=TRUE)],EstMGHsym$lambda,EstMGHsym$eta,0,0),4)
t3[,4] <- round(c(EstMGHsym$ep,NA,NA),4)

#NIG
t3[,5] <- round(c(c(EstNIG$beta1),c(EstNIG$sigmae),EstNIG$D1[upper.tri(EstNIG$D1,diag=TRUE)],EstNIG$lambda,EstNIG$eta,c(EstNIG$Gama)),4)
t3[,6] <- round(c(EstNIG$ep[1:7],NA,EstNIG$ep[8:10]),4)

t3

#################
## Table 4 
#################
t4 <- matrix(NA,ncol=5,nrow=3)
rownames(t4) <- c("loglik","AIC","BIC")
colnames(t4) <- c("GH", "GH sym", "NIG", "SN", "ST")

#MGH
t4[1,1] <- round(EstMGH$loglik,4)
t4[2,1] <- round(EstMGH$AIC,4)
t4[3,1] <- round(EstMGH$BIC,4)

#MGH sym
t4[1,2] <- round(EstMGHsym$loglik,4)
t4[2,2] <- round(EstMGHsym$AIC,4)
t4[3,2] <- round(EstMGHsym$BIC,4)

#MGH sym
t4[1,3] <- round(EstNIG$loglik,4)
t4[2,3] <- round(EstNIG$AIC,4)
t4[3,3] <- round(EstNIG$BIC,4)

#SN
t4[1,4] <- round(EstSn$loglik,4)
t4[2,4] <- round(EstSn$criteria$AIC,4)
t4[3,4] <- round(EstSn$criteria$BIC,4)

#ST
t4[1,5] <- round(EstST$loglik,4)
t4[2,5] <- round(EstST$criteria$AIC,4)
t4[3,5] <- round(EstST$criteria$BIC,4)

t4

#################################################################################################################################
## Plot residuals envelope
#################################################################################################################################
#pacotes
library(mvtnorm)
library(dplyr)
library(scales)

envelope(EstSn,nj,y,x,z,rep=500,level=.99)
envelope(EstST,nj,y,x,z,rep=500,level=.99)
envelope(EstMGH,nj,y,x,z,rep=500,level=.99)
envelope(EstNIG,nj,y,x,z,rep=500,level=.99)
