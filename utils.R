################################################################################
## likelihood
################################################################################

logvero<-function(nj,y,x,z,beta1,sigmae, D1, Gama, lambda, eta, type="MGH"){
  
  m<-length(nj)
  N<-sum(nj)
  p<-dim(x)[2]
  q1<-dim(z)[2]
  
  
  suma1<-0
  
  
  if (type=="MGH"){
    
    for (j in 1:m){
      y1<-y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))  ]
      x1<-matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
      z1<-matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
      med<-x1%*%beta1
      GamaF<-z1%*%Gama
      Psi<-(z1)%*%(D1)%*%t(z1)+diag(nj[j])*sigmae
      para<- ghyp(lambda = lambda, chi = eta, psi = eta, mu = med, sigma = Psi,
                  gamma = GamaF)
      suma1<-suma1+dghyp(y1, object = para, logvalue = TRUE)
    }
  }
  
  return(suma1)
  
}

################################################################################
## Funcao para maximizar
################################################################################

FCi<-function(pis, vj,ivj,lvj){
  m<-length(vj)
  lambda<-pis[1]
  eta<-pis[2]
  soma<-(lambda-1)*sum(lvj)-m*log(besselK(eta, lambda))-0.5*eta*sum(vj+ivj)
  return(-soma)
}

FCilambdafix<-function(eta,lambda, vj,ivj,lvj){
  m<-length(vj)
  #lambda<-pis[1]
  #eta<-pis[2]
  soma<-(lambda-1)*sum(lvj)-m*log(besselK(eta, lambda))-0.5*eta*sum(vj+ivj)
  return(-soma)
}

FCinig<-function(eta, vj,ivj,lvj){
  m<-length(vj)
  soma<-(-1.5)*sum(lvj)-m*log(besselK(eta, -0.5))-0.5*eta*sum(vj+ivj)
  return(-soma)
}

FCst<-function(eta, vj,ivj,lvj){
  m<-length(vj)
  soma<-(-eta/2-1)*sum(lvj)-m*log(besselK(eta, -eta/2))-0.5*eta*sum(vj+ivj)
  return(-soma)
}
