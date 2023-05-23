source("utils.R")
################################################################################
## The EM algorithm for Multivariate Generalized Hyperbolic LMM 
################################################################################
  
EM.MGH<-function(nj,y,x,z,beta1,sigmae,D1, Gama, lambda, eta,  type="MGH", precision=0.0000001, MaxIter=50){

###Instructions#########################################################################################################################################
# MLE of theta: beta(fixed effects), sigmae(error variance), D1(covariance matrix of random effects), lambda(skewness vector) and nu(vector or scalar).#
# nj: vector with size of clusters (vector)                                                                                                            #
# y: vector of responses (size = sum(nj) x 1)                                                                                                          #
# x: Design matrix of fixed effects (size = sum(nj) x p)                                                                                               #
# z: Design matrix of random effects (size = sum(nj) x q)                                                                                              #                                                                                                                                                              # 
#                                                                                                                         #
# precision: stopping rule for the maximisation                                                                                                        #                                                                                                     #
########################################################################################################################################################

	m<-length(nj)
	N<-sum(nj)
	p<-dim(x)[2]
	q1<-dim(z)[2]
	
	if (type=="MGH"){

	loglik<-1
	criterio<-1
	count<-0

	phis<-c(lambda,eta)

          
 teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],Gama, lambda,eta)

	while(criterio > precision){
		vj <- matrix(0,nrow=m,ncol=1)  # power 1
		ivj <- matrix(0,nrow=m,ncol=1)  # power -1
		lvj <- matrix(0,nrow=m,ncol=1)   # log vi
   
    sum1<- matrix(0,nrow=p,ncol=p)
		sum2<- matrix(0,nrow=p,ncol=1)
		sum3<-0
		sum4<-matrix(0,nrow=q1,ncol=q1)
		sum5<-matrix(0,nrow=q1,ncol=1)
    sum6<-matrix(0,nrow=q1,ncol=q1)
		bivj<-matrix(0,nrow=q1,ncol=m)
		bi<-matrix(0,nrow=q1,ncol=m)
   
  	Omega<-matrix(0,nrow=m*q1,ncol=q1)
		count <- count + 1
		cat("Iteration ",count,"\r")
		
    for (j in 1:m ){
		
    	y1<-y[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ]
			x1<-matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
			z1<-matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
			med<-x1%*%beta1
      GamaF<-z1%*%Gama
			Psi<-(z1)%*%(D1)%*%t(z1)+diag(nj[j])*sigmae
      Psi<-(Psi+t(Psi))/2
      iPsi<-solve(Psi)
			dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
      LambdaA<-solve(solve(D1)+t(z1)%*%z1/sigmae)
			Ai<-(diag(q1)-LambdaA%*%t(z1)%*%z1/sigmae)%*%Gama
      Bi<-LambdaA%*%t(z1)%*%(y1-med)/sigmae
      DeltaF<-eta+t(GamaF)%*%iPsi%*%GamaF
      
			vj[j] <- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("x"), check.pars = TRUE)
	    ivj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("1/x"), check.pars = TRUE)
    	lvj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("logx"), check.pars = TRUE)
    	##fix for NaN is expected value of logx
    	if(is.nan(lvj[j])){
    	   aux <- rgig(n = 10000, lambda = lambda-nj[j]/2 , chi = dj+eta, psi = DeltaF)
    	   lvj[j] <- mean(log(aux),na.rm=TRUE)
    	}
		
      bi[,j]<- Ai*vj[j]+Bi
      bivj[,j]<- Ai+Bi*ivj[j]      		

      Omega[((j-1)*q1+1):(j*q1),]<- vj[j]*(Ai%*%t(Ai))+(Ai%*%t(Bi)+Bi%*%t(Ai)+LambdaA)+ ivj[j]*(Bi%*%t(Bi))
		
			sum1<-sum1+(t(x1)%*%x1)*ivj[j]
			sum2<-sum2+(t(x1)%*%(y1*ivj[j]-z1%*%bivj[,j]))
      sum3<-sum3+ivj[j]*(t(y1-med)%*%(y1-med))-t(y1-med)%*%z1%*%bivj[,j]-t(bivj[,j])%*%t(z1)%*%(y1-med)+sum(diag(Omega[((j-1)*q1+1):(j*q1),]%*%t(z1)%*%z1))
      
      
			sum4<-sum4+Omega[((j-1)*q1+1):(j*q1),]-bi[,j]%*%t(Gama)-Gama%*%t(bi[,j])+vj[j]*Gama%*%t(Gama)
			sum5<- sum5 + solve(D1)%*%bi[,j]
      sum6<-sum6 + solve(D1)*vj[j]
		}

		beta1<-solve(sum1)%*%sum2
		sigmae<-as.numeric(sum3)/N
		Gama<-solve(sum6)%*%sum5
		D1<- sum4/m
		
    phis <- optim(phis, method = "L-BFGS-B", FCi, lower =c(-10,0.01), upper =c(10,10), vj=vj,ivj=ivj,lvj=lvj,hessian=TRUE)$par

    lambda<-phis[1]
    eta<-phis[2]
 
    dd<-D1[upper.tri(D1, diag = T)]

		lista<-list("beta"=beta1,"sigmae"=sigmae,"D"=dd,"Gama"=Gama,"lambda"=lambda,"eta"=eta,"iter"=count,"criterio"=criterio)
		
		logvero1<-logvero(nj,y,x,z,beta1,sigmae, D1, Gama, lambda, eta, type="MGH")

		loglik1 <- loglik
		loglik<- logvero1

		criterio <- abs(loglik/loglik1-1)

		if (count==MaxIter){
			criterio <- 1e-16
		}
	}
  teta <- c(beta1, sigmae, D1[upper.tri(D1, diag = T)],Gama, eta)
 }  ## FiM IF

 	if (type=="MGHsym"){ ### gamma=0

	loglik<-1
	criterio<-1
	count<-0

	phis<-c(lambda,eta)

	while(criterio > precision){
		vj <- NULL   # power 1
		ivj <- NULL  # power -1
		lvj <- NULL   # log vi
   
    sum1<- matrix(0,nrow=p,ncol=p)
		sum2<- matrix(0,nrow=p,ncol=1)
		sum3<-0
		sum4<-matrix(0,nrow=q1,ncol=q1)
		sum5<-matrix(0,nrow=q1,ncol=1)
		bivj<-matrix(0,nrow=q1,ncol=m)
		bi<-matrix(0,nrow=q1,ncol=m)
   
  	Omega<-matrix(0,nrow=m*q1,ncol=q1)
		count <- count + 1
		cat("Iteration ",count,"\r")
		
    for (j in 1:m ){
		
    	y1<-y[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ]
			x1<-matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
			z1<-matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
			med<-x1%*%beta1
      GamaF<-z1%*%Gama
			Psi<-(z1)%*%(D1)%*%t(z1)+diag(nj[j])*sigmae
      Psi<-(Psi+t(Psi))/2
      iPsi<-solve(Psi)
			dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
      LambdaA<-solve(solve(D1)+t(z1)%*%z1/sigmae)
			Ai<-(diag(q1)-LambdaA%*%t(z1)%*%z1/sigmae)%*%Gama
      Bi<-LambdaA%*%t(z1)%*%(y1-med)/sigmae
      DeltaF<-eta+t(GamaF)%*%iPsi%*%GamaF
      
			vj[j] <- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("x"), check.pars = TRUE)
	    ivj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("1/x"), check.pars = TRUE)
    	lvj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("logx"), check.pars = TRUE)
    	##fix for NaN is expected value of logx
    	if(is.nan(lvj[j])){
    	   aux <- rgig(n = 10000, lambda = lambda-nj[j]/2 , chi = dj+eta, psi = DeltaF)
    	   lvj[j] <- mean(log(aux),na.rm=TRUE)
    	}

			bi[,j]<- Ai*vj[j]+Bi
      bivj[,j]<- Ai+Bi*ivj[j]      		

      Omega[((j-1)*q1+1):(j*q1),]<- vj[j]*(Ai%*%t(Ai))+(Ai%*%t(Bi)+Bi%*%t(Ai)+LambdaA)+ ivj[j]*(Bi%*%t(Bi))
		
			sum1<-sum1+(t(x1)%*%x1)*ivj[j]
			sum2<-sum2+(t(x1)%*%(y1*ivj[j]-z1%*%bivj[,j]))
      sum3<-sum3+ivj[j]*(t(y1-med)%*%(y1-med))-t(y1-med)%*%z1%*%bivj[,j]-t(bivj[,j])%*%t(z1)%*%(y1-med)+sum(diag(Omega[((j-1)*q1+1):(j*q1),]%*%t(z1)%*%z1))
      
      
			sum4<-sum4+Omega[((j-1)*q1+1):(j*q1),]-bi[,j]%*%t(Gama)-Gama%*%t(bi[,j])+vj[j]*Gama%*%t(Gama)
			sum5<- sum5 + bi[,j]
		}

		beta1<-solve(sum1)%*%sum2
		sigmae<-as.numeric(sum3)/N
		Gama<-matrix(0,nrow=q1,ncol=1)
		D1<- sum4/m
		
		phis <- optim(phis, method = "L-BFGS-B", FCi, lower =c(-10,0.01), upper =c(10,200), vj=vj,ivj=ivj,lvj=lvj,hessian=TRUE)$par
    lambda<-phis[1]
    eta<-phis[2]
 
    dd<-D1[upper.tri(D1, diag = T)]

		lista<-list("beta"=beta1,"sigmae"=sigmae,"D"=dd,"Gama"=Gama,"lambda"=lambda,"eta"=eta,"iter"=count,"criterio"=criterio)
		
		logvero1<-logvero(nj,y,x,z,beta1,sigmae, D1, Gama, lambda, eta, type="MGH")

		loglik1 <- loglik
		loglik<- logvero1

		criterio <- abs(loglik/loglik1-1)

		if (count==MaxIter){
			criterio <- 1e-16
		}
	}
  teta <- c(beta1, sigmae, D1[upper.tri(D1, diag = T)], lambda, eta)
 }  ## FiM IF
 
 	if (type=="NIG"){ ### lambda= -0.5

	loglik<-1
	criterio<-1
	count<-0

	while(criterio > precision){
		vj <- NULL   # power 1
		ivj <- NULL  # power -1
		lvj <- NULL   # log vi
   
    sum1<- matrix(0,nrow=p,ncol=p)
		sum2<- matrix(0,nrow=p,ncol=1)
		sum3<-0
		sum4<-matrix(0,nrow=q1,ncol=q1)
		sum5<-matrix(0,nrow=q1,ncol=1)
		bivj<-matrix(0,nrow=q1,ncol=m)
		bi<-matrix(0,nrow=q1,ncol=m)
   
  	Omega<-matrix(0,nrow=m*q1,ncol=q1)
		count <- count + 1
		cat("Iteration ",count,"\r")
		
    for (j in 1:m ){
		
    	y1<-y[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ]
			x1<-matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
			z1<-matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
			med<-x1%*%beta1
      GamaF<-z1%*%Gama
			Psi<-(z1)%*%(D1)%*%t(z1)+diag(nj[j])*sigmae
      Psi<-(Psi+t(Psi))/2
      iPsi<-solve(Psi)
			dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
      LambdaA<-solve(solve(D1)+t(z1)%*%z1/sigmae)
			Ai<-(diag(q1)-LambdaA%*%t(z1)%*%z1/sigmae)%*%Gama
      Bi<-LambdaA%*%t(z1)%*%(y1-med)/sigmae
      DeltaF<-eta+t(GamaF)%*%iPsi%*%GamaF
      
			vj[j] <- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("x"), check.pars = TRUE)
	    ivj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("1/x"), check.pars = TRUE)
    	lvj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("logx"), check.pars = TRUE)
      			
			bi[,j]<- Ai*vj[j]+Bi
      bivj[,j]<- Ai+Bi*ivj[j]      		

      Omega[((j-1)*q1+1):(j*q1),]<- vj[j]*(Ai%*%t(Ai))+(Ai%*%t(Bi)+Bi%*%t(Ai)+LambdaA)+ ivj[j]*(Bi%*%t(Bi))
		
			sum1<-sum1+(t(x1)%*%x1)*ivj[j]
			sum2<-sum2+(t(x1)%*%(y1*ivj[j]-z1%*%bivj[,j]))
      sum3<-sum3+ivj[j]*(t(y1-med)%*%(y1-med))-t(y1-med)%*%z1%*%bivj[,j]-t(bivj[,j])%*%t(z1)%*%(y1-med)+sum(diag(Omega[((j-1)*q1+1):(j*q1),]%*%t(z1)%*%z1))
      
      
			sum4<-sum4+Omega[((j-1)*q1+1):(j*q1),]-bi[,j]%*%t(Gama)-Gama%*%t(bi[,j])+vj[j]*Gama%*%t(Gama)
			sum5<- sum5 + bi[,j]
		}

		beta1<-solve(sum1)%*%sum2
		sigmae<-as.numeric(sum3)/N
		Gama<-sum5/as.numeric(sum(vj))
		D1<- sum4/m
		
		eta <- optim(eta, method = "L-BFGS-B", FCinig, lower =0.01, upper =200, vj=vj,ivj=ivj,lvj=lvj,hessian=TRUE)$par
    lambda<- -0.5
    
    dd<-D1[upper.tri(D1, diag = T)]

		lista<-list("beta"=beta1,"sigmae"=sigmae,"D"=dd,"Gama"=Gama,"lambda"=lambda,"eta"=eta,"iter"=count,"criterio"=criterio)
		
		logvero1<-logvero(nj,y,x,z,beta1,sigmae, D1, Gama, lambda, eta, type="MGH")

		loglik1 <- loglik
		loglik<- logvero1

		criterio <- abs(loglik/loglik1-1)

		if (count==MaxIter){
			criterio <- 1e-16
		}
	}
  teta <- c(beta1, sigmae, D1[upper.tri(D1, diag = T)],Gama,  eta)
 }  ## FiM IF

	npar<-length(teta)

	##Model comparison criteria
	AICc<- -2*loglik +2*npar
	AICcorr<- AICc + ((2*npar*(npar+1))/(sum(nj)-npar-1))
	BICc <- -2*loglik +log(sum(nj))*npar
 

 ep<-matrizMI(nj, y, x, z, beta1, sigmae, D1, Gama, lambda, eta , type=type)
 beta_nam <- paste0("beta_",0:(length(beta1)-1))
 D_nam <- paste0("D_1",1:ncol(D1))
 for(i in 2:nrow(D1)) D_nam <- c(D_nam,paste0("D_",i,i:ncol(D1)))
 names <- c(beta_nam,"sigma_e",D_nam)
 if(type=="MGH") names <- c(names,"gamma_1","gamma_2","lambda","eta")
 if(type=="MGHsym") names <- c(names,"lambda","eta")
 if(type=="NIG") names <- c(names,"gamma_1","gamma_2","eta")
 names(ep) <- names
 obj.out <- list(beta1=beta1,sigmae=sigmae, D1=D1, lambda=lambda, eta = eta, Gama = Gama, type=type,  loglik=loglik, AIC=AICc, AICc=AICcorr, BIC=BICc,  iter = count,ep=ep)
	

	class(obj.out) <- "EM.MGH"
	return(obj.out)
}


################################################################################
## The EM algorithm for Multivariate Generalized Hyperbolic LMM 
################################################################################

EM.MGHini<-function(nj,y,x,z,beta1,sigmae,D1, Gama, lambda, eta,  type="MGH",MaxIter=50){
  
  ###Instructions#########################################################################################################################################
  # MLE of theta: beta(fixed effects), sigmae(error variance), D1(covariance matrix of random effects), lambda(skewness vector) and nu(vector or scalar).#
  # nj: vector with size of clusters (vector)                                                                                                            #
  # y: vector of responses (size = sum(nj) x 1)                                                                                                          #
  # x: Design matrix of fixed effects (size = sum(nj) x p)                                                                                               #
  # z: Design matrix of random effects (size = sum(nj) x q)                                                                                              #                                                                                                                                                              # 
  #                                                                                                                         #
  # precision: stopping rule for the maximisation                                                                                                        #                                                                                                     #
  ########################################################################################################################################################
  
  m<-length(nj)
  N<-sum(nj)
  p<-dim(x)[2]
  q1<-dim(z)[2]
  
  if (type=="MGH"){
    
    loglik<-1
    count<-0
    
    phis<-c(lambda,eta)
    
    teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],Gama, lambda,eta)
    
    while(count < MaxIter){
      vj <- matrix(0,nrow=m,ncol=1)  # power 1
      ivj <- matrix(0,nrow=m,ncol=1)  # power -1
      lvj <- matrix(0,nrow=m,ncol=1)   # log vi
      
      sum1<- matrix(0,nrow=p,ncol=p)
      sum2<- matrix(0,nrow=p,ncol=1)
      sum3<-0
      sum4<-matrix(0,nrow=q1,ncol=q1)
      sum5<-matrix(0,nrow=q1,ncol=1)
      sum6<-matrix(0,nrow=q1,ncol=q1)
      bivj<-matrix(0,nrow=q1,ncol=m)
      bi<-matrix(0,nrow=q1,ncol=m)
      
      Omega<-matrix(0,nrow=m*q1,ncol=q1)
      count <- count + 1
      cat("Iteration ",count,"\r")
      
      for (j in 1:m ){
        
        y1<-y[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ]
        x1<-matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
        z1<-matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
        med<-x1%*%beta1
        GamaF<-z1%*%Gama
        Psi<-(z1)%*%(D1)%*%t(z1)+diag(nj[j])*sigmae
        Psi<-(Psi+t(Psi))/2
        iPsi<-solve(Psi)
        dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
        LambdaA<-solve(solve(D1)+t(z1)%*%z1/sigmae)
        Ai<-(diag(q1)-LambdaA%*%t(z1)%*%z1/sigmae)%*%Gama
        Bi<-LambdaA%*%t(z1)%*%(y1-med)/sigmae
        DeltaF<-eta+t(GamaF)%*%iPsi%*%GamaF
        
        vj[j] <- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("x"), check.pars = TRUE)
        ivj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("1/x"), check.pars = TRUE)
        lvj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("logx"), check.pars = TRUE)
        
        ##fix for NaN is expected value of logx
        if(is.nan(lvj[j])){
          aux <- rgig(n = 10000, lambda = lambda-nj[j]/2 , chi = dj+eta, psi = DeltaF)
          lvj[j] <- mean(log(aux),na.rm=TRUE)
        }
        
        bi[,j]<- Ai*vj[j]+Bi
        bivj[,j]<- Ai+Bi*ivj[j]      		
        
        Omega[((j-1)*q1+1):(j*q1),]<- vj[j]*(Ai%*%t(Ai))+(Ai%*%t(Bi)+Bi%*%t(Ai)+LambdaA)+ ivj[j]*(Bi%*%t(Bi))
        
        sum1<-sum1+(t(x1)%*%x1)*ivj[j]
        sum2<-sum2+(t(x1)%*%(y1*ivj[j]-z1%*%bivj[,j]))
        sum3<-sum3+ivj[j]*(t(y1-med)%*%(y1-med))-t(y1-med)%*%z1%*%bivj[,j]-t(bivj[,j])%*%t(z1)%*%(y1-med)+sum(diag(Omega[((j-1)*q1+1):(j*q1),]%*%t(z1)%*%z1))
        
        sum4<-sum4+Omega[((j-1)*q1+1):(j*q1),]-bi[,j]%*%t(Gama)-Gama%*%t(bi[,j])+vj[j]*Gama%*%t(Gama)
        sum5<- sum5 + solve(D1)%*%bi[,j]
        sum6<-sum6 + solve(D1)*vj[j]
      }
      
      beta1<-solve(sum1)%*%sum2
      sigmae<-as.numeric(sum3)/N
      Gama<-solve(sum6)%*%sum5
      D1<- sum4/m
      
      phis <- optim(phis, method = "L-BFGS-B", FCi, lower =c(-100,0.01), upper =c(100,100), vj=vj,ivj=ivj,lvj=lvj,hessian=TRUE)$par
      
      lambda<-phis[1]
      eta<-phis[2]
      
      dd<-D1[upper.tri(D1, diag = T)]
      
      lista<-list("beta"=beta1,"sigmae"=sigmae,"D"=dd,"Gama"=Gama,"lambda"=lambda,"eta"=eta,"iter"=count)
      
      logvero1<-logvero(nj,y,x,z,beta1,sigmae, D1, Gama, lambda, eta, type="MGH")
      
      loglik1 <- loglik
      loglik<- logvero1
      
      param <- teta
    
    }
    teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],Gama,lambda,eta)
  }  ## FiM IF
  
  if (type=="MGHsym"){ ### gamma=0
    
    loglik<-1
    count<-0
    
    phis<-c(lambda,eta)
    
    while(count < MaxIter){
      vj <- NULL   # power 1
      ivj <- NULL  # power -1
      lvj <- NULL   # log vi
      
      sum1<- matrix(0,nrow=p,ncol=p)
      sum2<- matrix(0,nrow=p,ncol=1)
      sum3<-0
      sum4<-matrix(0,nrow=q1,ncol=q1)
      sum5<-matrix(0,nrow=q1,ncol=1)
      bivj<-matrix(0,nrow=q1,ncol=m)
      bi<-matrix(0,nrow=q1,ncol=m)
      
      Omega<-matrix(0,nrow=m*q1,ncol=q1)
      count <- count + 1
      cat("Iteration ",count,"\r")
      
      for (j in 1:m ){
        
        y1<-y[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ]
        x1<-matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
        z1<-matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
        med<-x1%*%beta1
        GamaF<-z1%*%Gama
        Psi<-(z1)%*%(D1)%*%t(z1)+diag(nj[j])*sigmae
        Psi<-(Psi+t(Psi))/2
        iPsi<-solve(Psi)
        dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
        LambdaA<-solve(solve(D1)+t(z1)%*%z1/sigmae)
        Ai<-(diag(q1)-LambdaA%*%t(z1)%*%z1/sigmae)%*%Gama
        Bi<-LambdaA%*%t(z1)%*%(y1-med)/sigmae
        DeltaF<-eta+t(GamaF)%*%iPsi%*%GamaF
        
        vj[j] <- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("x"), check.pars = TRUE)
        ivj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("1/x"), check.pars = TRUE)
        lvj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("logx"), check.pars = TRUE)
        
        bi[,j]<- Ai*vj[j]+Bi
        bivj[,j]<- Ai+Bi*ivj[j]      		
        
        Omega[((j-1)*q1+1):(j*q1),]<- vj[j]*(Ai%*%t(Ai))+(Ai%*%t(Bi)+Bi%*%t(Ai)+LambdaA)+ ivj[j]*(Bi%*%t(Bi))
        
        sum1<-sum1+(t(x1)%*%x1)*ivj[j]
        sum2<-sum2+(t(x1)%*%(y1*ivj[j]-z1%*%bivj[,j]))
        sum3<-sum3+ivj[j]*(t(y1-med)%*%(y1-med))-t(y1-med)%*%z1%*%bivj[,j]-t(bivj[,j])%*%t(z1)%*%(y1-med)+sum(diag(Omega[((j-1)*q1+1):(j*q1),]%*%t(z1)%*%z1))
        
        
        sum4<-sum4+Omega[((j-1)*q1+1):(j*q1),]-bi[,j]%*%t(Gama)-Gama%*%t(bi[,j])+vj[j]*Gama%*%t(Gama)
        sum5<- sum5 + bi[,j]
      }
      
      beta1<-solve(sum1)%*%sum2
      sigmae<-as.numeric(sum3)/N
      Gama<-matrix(0,nrow=q1,ncol=1)
      D1<- sum4/m
      
      phis <- optim(phis, method = "L-BFGS-B", FCi, lower =c(-10,0.01), upper =c(10,200), vj=vj,ivj=ivj,lvj=lvj,hessian=TRUE)$par
      lambda<-phis[1]
      eta<-phis[2]
      
      dd<-D1[upper.tri(D1, diag = T)]
      
      lista<-list("beta"=beta1,"sigmae"=sigmae,"D"=dd,"Gama"=Gama,"lambda"=lambda,"eta"=eta,"iter"=count)
      
      logvero1<-logvero(nj,y,x,z,beta1,sigmae, D1, Gama, lambda, eta, type="MGH")
      
      loglik1 <- loglik
      loglik<- logvero1
      
    }
    teta <- c(beta1, sigmae, D1[upper.tri(D1, diag = T)], lambda, eta)
  }  ## FiM IF
  
  if (type=="NIG"){ ### lambda= -0.5
    
    loglik<-1
    count<-0
    
    while(count < MaxIter){
      vj <- NULL   # power 1
      ivj <- NULL  # power -1
      lvj <- NULL   # log vi
      
      sum1<- matrix(0,nrow=p,ncol=p)
      sum2<- matrix(0,nrow=p,ncol=1)
      sum3<-0
      sum4<-matrix(0,nrow=q1,ncol=q1)
      sum5<-matrix(0,nrow=q1,ncol=1)
      bivj<-matrix(0,nrow=q1,ncol=m)
      bi<-matrix(0,nrow=q1,ncol=m)
      
      Omega<-matrix(0,nrow=m*q1,ncol=q1)
      count <- count + 1
      cat("Iteration ",count,"\r")
      
      for (j in 1:m ){
        
        y1<-y[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ]
        x1<-matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
        z1<-matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
        med<-x1%*%beta1
        GamaF<-z1%*%Gama
        Psi<-(z1)%*%(D1)%*%t(z1)+diag(nj[j])*sigmae
        Psi<-(Psi+t(Psi))/2
        iPsi<-solve(Psi)
        dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
        LambdaA<-solve(solve(D1)+t(z1)%*%z1/sigmae)
        Ai<-(diag(q1)-LambdaA%*%t(z1)%*%z1/sigmae)%*%Gama
        Bi<-LambdaA%*%t(z1)%*%(y1-med)/sigmae
        DeltaF<-eta+t(GamaF)%*%iPsi%*%GamaF
        
        vj[j] <- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("x"), check.pars = TRUE)
        ivj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("1/x"), check.pars = TRUE)
        lvj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("logx"), check.pars = TRUE)
        
        bi[,j]<- Ai*vj[j]+Bi
        bivj[,j]<- Ai+Bi*ivj[j]      		
        
        Omega[((j-1)*q1+1):(j*q1),]<- vj[j]*(Ai%*%t(Ai))+(Ai%*%t(Bi)+Bi%*%t(Ai)+LambdaA)+ ivj[j]*(Bi%*%t(Bi))
        
        sum1<-sum1+(t(x1)%*%x1)*ivj[j]
        sum2<-sum2+(t(x1)%*%(y1*ivj[j]-z1%*%bivj[,j]))
        sum3<-sum3+ivj[j]*(t(y1-med)%*%(y1-med))-t(y1-med)%*%z1%*%bivj[,j]-t(bivj[,j])%*%t(z1)%*%(y1-med)+sum(diag(Omega[((j-1)*q1+1):(j*q1),]%*%t(z1)%*%z1))
        
        
        sum4<-sum4+Omega[((j-1)*q1+1):(j*q1),]-bi[,j]%*%t(Gama)-Gama%*%t(bi[,j])+vj[j]*Gama%*%t(Gama)
        sum5<- sum5 + bi[,j]
      }
      
      beta1<-solve(sum1)%*%sum2
      sigmae<-as.numeric(sum3)/N
      Gama<-sum5/as.numeric(sum(vj))
      D1<- sum4/m
      
      eta <- optim(eta, method = "L-BFGS-B", FCinig, lower =0.01, upper =200, vj=vj,ivj=ivj,lvj=lvj,hessian=TRUE)$par
      lambda<- -0.5
      
      dd<-D1[upper.tri(D1, diag = T)]
      
      lista<-list("beta"=beta1,"sigmae"=sigmae,"D"=dd,"Gama"=Gama,"lambda"=lambda,"eta"=eta,"iter"=count)
      
      logvero1<-logvero(nj,y,x,z,beta1,sigmae, D1, Gama, lambda, eta, type="MGH")
      
      
      loglik1 <- loglik
      loglik<- logvero1
    }
    teta <- c(beta1, sigmae, D1[upper.tri(D1, diag = T)],Gama,  eta)
  }  ## FiM IF
  
  obj.out <- list(beta1=beta1,sigmae=sigmae, D1=D1, lambda=lambda, eta = eta, Gama = Gama, type=type,  loglik=loglik, iter = count)
  
  class(obj.out) <- "EM.MGH.INI"
  return(obj.out)
}



 
