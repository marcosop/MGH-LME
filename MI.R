derBK<-function(x,lambda){
 h <- 1e-10
ff <- function(x,lambda){besselK(x, lambda)}
aa<-(ff(x,lambda+h) - ff(x,lambda)) / h
return(aa)
} 




################################################################
#Information Matrix
################################################################
      
matrizMI<-function(nj, y, x, z, beta1, sigmae, D1, Gama, lambda, eta , type="MGH") {

# MLE of theta: beta1(vector),sigmae(random errors),D1(matrix),Gama (vector),lambda(scalar) and eta(scalar).
# nj: vector with size of clusters (vetor )
# y: vector of Responses      (size=sum(nj) x 1)
# x: Design matrix of fixed effects      (size=sum(nj) x p)
# z: Design matrix of random effects      (size=sum(nj) x q)
m<-length(nj)
N<-sum(nj)
p<-dim(x)[2]
q1<-dim(z)[2]

    vj<-  matrix(0,nrow=m,ncol=1)   # power 1
		ivj<- matrix(0,nrow=m,ncol=1)  # power -1
		lvj<- matrix(0,nrow=m,ncol=1)   # log vi
   
   	bivj<-matrix(0,nrow=q1,ncol=m)
		bi<-matrix(0,nrow=q1,ncol=m)
   
  	Omega<-matrix(0,nrow=m*q1,ncol=q1)

if (type=="MGH"){

npar<-p+1+q1*(q1+1)/2+q1+2

suma1<-matrix(0,npar,npar)

for (j in 1:m){        
      y1<- y[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ]
			x1<- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
			z1<- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
			med<- x1%*%beta1
      GamaF<-z1%*%Gama
			Psi<-(z1)%*%(D1)%*%t(z1)+diag(nj[j])*sigmae
      Psi<-(Psi+t(Psi))/2
      iPsi<-solve(Psi,tol=1e-30)
			dj<-as.numeric(t(y1-med)%*%iPsi%*%(y1-med))
      LambdaA<-solve(solve(D1,tol=1e-30)+t(z1)%*%z1/sigmae,tol=1e-30)
			Ai<-(diag(q1)-LambdaA%*%t(z1)%*%z1/sigmae)%*%Gama
      Bi<-LambdaA%*%t(z1)%*%(y1-med)/sigmae
      DeltaF<-as.numeric(eta+t(GamaF)%*%iPsi%*%GamaF)
      
      omegai<-as.numeric(sqrt((eta+dj)*(DeltaF)))
      omegaibeta<- -as.numeric(omegai*DeltaF)*(t(x1)%*%iPsi%*%(y1-med))
			
      vj[j] <- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("x"))
	    ivj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("1/x"))
    	lvj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("logx"))
    	if(is.nan(lvj[j])){
    	   #print(c(lvj[j],lambda, dj, eta, DeltaF))
    	   aux <- rgig(n = 10000, lambda = lambda-nj[j]/2 , chi = dj+eta, psi = DeltaF)
    	   #print(aux)
    	   lvj[j] <- mean(log(aux),na.rm=TRUE)
    	}

      bi[,j]<- Ai*vj[j]+Bi
      bivj[,j]<- Ai+Bi*ivj[j]      		

      Omega[((j-1)*q1+1):(j*q1),]<- vj[j]*(Ai%*%t(Ai))+(Ai%*%t(Bi)+Bi%*%t(Ai)+LambdaA)+ ivj[j]*(Bi%*%t(Bi))
      
        derAjbeta<- ( (t(x1)%*%(y1-med))*ivj[j]-t(x1)%*%z1%*%bivj[,j] )/sigmae
        derAjsigmae <- (ivj[j]*(t(y1-med)%*%(y1-med))-t(y1-med)%*%z1%*%bivj[,j]-t(bivj[,j])%*%t(z1)%*%(y1-med)+sum(diag(Omega[((j-1)*q1+1):(j*q1),]%*%t(z1)%*%z1)))/(2*sigmae^2)-nj[j]/(2*sigmae)
        
        Baa<-matrix(c(1, 0, 0 ,0),ncol=2,nrow=2)
        Caa<-matrix(c(0, 1,1,0),nrow=2,ncol=2)
        Daa<-matrix(c(0,0,0,1),nrow=2,ncol=2)
        Aux1<-Omega[((j-1)*q1+1):(j*q1),]-bi[,j]%*%t(Gama)-Gama%*%t(bi[,j])+vj[j]*Gama%*%t(Gama)
        derAjalfa1<- -0.5*sum(diag(solve(D1)%*%Baa))+0.5*sum(diag(solve(D1)%*%Baa%*%solve(D1)%*%(Aux1)))
        derAjalfa2<- -0.5*sum(diag(solve(D1)%*%Caa))+0.5*sum(diag(solve(D1)%*%Caa%*%solve(D1)%*%(Aux1)))
        derAjalfa3<- -0.5*sum(diag(solve(D1)%*%Daa))+0.5*sum(diag(solve(D1)%*%Daa%*%solve(D1)%*%(Aux1)))

        derAjgama<- solve(D1)%*%bi[,j]-vj[j]*(solve(D1)%*%Gama)
        derAjlambda<- lvj[j]-derBK(eta,lambda)/besselK(eta, lambda)
        derAjeta<- (lambda/eta+besselK(eta, lambda-1)/besselK(eta, lambda))-0.5*(vj[j]+ivj[j])
        
        derlog<-  matrix(c(derAjbeta,derAjsigmae,derAjalfa1,derAjalfa2,derAjalfa3,derAjgama,derAjlambda,derAjeta),npar,1)
        suma1<- suma1+derlog%*%t(derlog)
}
  }
  
  
if (type=="MGHsym"){

npar<-p+1+q1*(q1+1)/2+2

suma1<-matrix(0,npar,npar)

for (j in 1:m){        
      y1<- y[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ]
			x1<- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
			z1<- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
			med<- x1%*%beta1
      GamaF<-z1%*%Gama
			Psi<-(z1)%*%(D1)%*%t(z1)+diag(nj[j])*sigmae
      Psi<-(Psi+t(Psi))/2
      iPsi<-solve(Psi,tol=1e-30)
			dj<-as.numeric(t(y1-med)%*%iPsi%*%(y1-med))
      LambdaA<-solve(solve(D1,tol=1e-30)+t(z1)%*%z1/sigmae,tol=1e-30)
			Ai<-(diag(q1)-LambdaA%*%t(z1)%*%z1/sigmae)%*%Gama
      Bi<-LambdaA%*%t(z1)%*%(y1-med)/sigmae
      DeltaF<-as.numeric(eta+t(GamaF)%*%iPsi%*%GamaF)
      
      omegai<-as.numeric(sqrt((eta+dj)*(DeltaF)))
      omegaibeta<- -as.numeric(omegai*DeltaF)*(t(x1)%*%iPsi%*%(y1-med))
			
      vj[j] <- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("x"))
      ivj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("1/x"))
      lvj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("logx"))
      ##fix for NaN is expected value of logx
    	if(is.nan(lvj[j])){
    	   #print(c(count, j, lambda, dj, eta, DeltaF))
    	   #print(c(lvj[j],lambda, dj, eta, DeltaF))
    	   aux <- rgig(n = 10000, lambda = lambda-nj[j]/2 , chi = dj+eta, psi = DeltaF)
    	   #print(aux)
    	   lvj[j] <- mean(log(aux),na.rm=TRUE)
    	}
			
			bi[,j]<- Ai*vj[j]+Bi
      bivj[,j]<- Ai+Bi*ivj[j]      		

      Omega[((j-1)*q1+1):(j*q1),]<- vj[j]*(Ai%*%t(Ai))+(Ai%*%t(Bi)+Bi%*%t(Ai)+LambdaA)+ ivj[j]*(Bi%*%t(Bi))
      
        derAjbeta<- ( (t(x1)%*%(y1-med))*ivj[j]-t(x1)%*%z1%*%bivj[,j] )/sigmae
        derAjsigmae <- (ivj[j]*(t(y1-med)%*%(y1-med))-t(y1-med)%*%z1%*%bivj[,j]-t(bivj[,j])%*%t(z1)%*%(y1-med)+sum(diag(Omega[((j-1)*q1+1):(j*q1),]%*%t(z1)%*%z1)))/(2*sigmae^2)-nj[j]/(2*sigmae)
        
        Baa<-matrix(c(1, 0, 0 ,0),ncol=2,nrow=2)
        Caa<-matrix(c(0, 1,1,0),nrow=2,ncol=2)
        Daa<-matrix(c(0,0,0,1),nrow=2,ncol=2)
        Aux1<-Omega[((j-1)*q1+1):(j*q1),]-bi[,j]%*%t(Gama)-Gama%*%t(bi[,j])+vj[j]*Gama%*%t(Gama)
        derAjalfa1<- -0.5*sum(diag(solve(D1)%*%Baa))+0.5*sum(diag(solve(D1)%*%Baa%*%solve(D1)%*%(Aux1)))
        derAjalfa2<- -0.5*sum(diag(solve(D1)%*%Caa))+0.5*sum(diag(solve(D1)%*%Caa%*%solve(D1)%*%(Aux1)))
        derAjalfa3<- -0.5*sum(diag(solve(D1)%*%Daa))+0.5*sum(diag(solve(D1)%*%Daa%*%solve(D1)%*%(Aux1)))

        #derAjgama<- solve(D1)%*%bi[,j]-vj[j]*(solve(D1)%*%Gama)
        derAjlambda<- lvj[j]-derBK(eta,lambda)/besselK(eta, lambda)
        derAjeta<- (lambda/eta+besselK(eta, lambda-1)/besselK(eta, lambda))-0.5*(vj[j]+ivj[j])
        
        derlog<-  matrix(c(derAjbeta,derAjsigmae,derAjalfa1,derAjalfa2,derAjalfa3,derAjlambda,derAjeta),npar,1)
        suma1<- suma1+derlog%*%t(derlog)
}
  }
  
 if (type=="NIG"){

npar<-p+1+q1*(q1+1)/2+q1+1

suma1<-matrix(0,npar,npar)

for (j in 1:m){        
      y1<- y[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ]
			x1<- matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
			z1<- matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
			med<- x1%*%beta1
      GamaF<-z1%*%Gama
			Psi<-(z1)%*%(D1)%*%t(z1)+diag(nj[j])*sigmae
      Psi<-(Psi+t(Psi))/2
      iPsi<-solve(Psi,tol=1e-30)
			dj<-as.numeric(t(y1-med)%*%iPsi%*%(y1-med))
      LambdaA<-solve(solve(D1,tol=1e-30)+t(z1)%*%z1/sigmae,tol=1e-30)
			Ai<-(diag(q1)-LambdaA%*%t(z1)%*%z1/sigmae)%*%Gama
      Bi<-LambdaA%*%t(z1)%*%(y1-med)/sigmae
      DeltaF<-as.numeric(eta+t(GamaF)%*%iPsi%*%GamaF)
      
    
      vj[j] <- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("x"))
	    ivj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("1/x"))
    	lvj[j]<- Egig(lambda-nj[j]/2, dj+eta, DeltaF, func = c("logx"))
      			
			bi[,j]<- Ai*vj[j]+Bi
      bivj[,j]<- Ai+Bi*ivj[j]      		

      Omega[((j-1)*q1+1):(j*q1),]<- vj[j]*(Ai%*%t(Ai))+(Ai%*%t(Bi)+Bi%*%t(Ai)+LambdaA)+ ivj[j]*(Bi%*%t(Bi))
      
        derAjbeta<- ( (t(x1)%*%(y1-med))*ivj[j]-t(x1)%*%z1%*%bivj[,j] )/sigmae
        derAjsigmae <- (ivj[j]*(t(y1-med)%*%(y1-med))-t(y1-med)%*%z1%*%bivj[,j]-t(bivj[,j])%*%t(z1)%*%(y1-med)+sum(diag(Omega[((j-1)*q1+1):(j*q1),]%*%t(z1)%*%z1)))/(2*sigmae^2)-nj[j]/(2*sigmae)
        
        Baa<-matrix(c(1, 0, 0 ,0),ncol=2,nrow=2)
        Caa<-matrix(c(0, 1,1,0),nrow=2,ncol=2)
        Daa<-matrix(c(0,0,0,1),nrow=2,ncol=2)
        Aux1<-Omega[((j-1)*q1+1):(j*q1),]-bi[,j]%*%t(Gama)-Gama%*%t(bi[,j])+vj[j]*Gama%*%t(Gama)
        derAjalfa1<- -0.5*sum(diag(solve(D1)%*%Baa))+0.5*sum(diag(solve(D1)%*%Baa%*%solve(D1)%*%(Aux1)))
        derAjalfa2<- -0.5*sum(diag(solve(D1)%*%Caa))+0.5*sum(diag(solve(D1)%*%Caa%*%solve(D1)%*%(Aux1)))
        derAjalfa3<- -0.5*sum(diag(solve(D1)%*%Daa))+0.5*sum(diag(solve(D1)%*%Daa%*%solve(D1)%*%(Aux1)))

        derAjgama<- solve(D1)%*%bi[,j]-vj[j]*(solve(D1)%*%Gama)
       # derAjlambda<- lvj[j]-derBK(eta,lambda)/besselK(eta, lambda)
        derAjeta<- (lambda/eta+besselK(eta, lambda-1)/besselK(eta, lambda))-0.5*(vj[j]+ivj[j])        
        derlog<-  matrix(c(derAjbeta,derAjsigmae,derAjalfa1,derAjalfa2,derAjalfa3,derAjgama,derAjeta),npar,1)
        suma1<- suma1+derlog%*%t(derlog)
}
  }
  
desvio<- sqrt(diag(solve(suma1)))
return(desvio)
#return(suma1)
}


#t(x1)%*%z1%*%bivj[,j]+t(bivj[,j])%*%t(z1)%*%(x1)
