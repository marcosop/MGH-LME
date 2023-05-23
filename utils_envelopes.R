logvero.lmm = function(y,x,z,ind,beta1,sigmae,D1,lambda,distr,nu){ #ind = indicadora de individuo
  m<-n_distinct(ind)
  N<-length(ind)
  p<-dim(x)[2]
  q1<-dim(z)[2]

  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda));
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
lv=c(0)

for (jj in 1:N){
  if (distr=="sn") lv[jj] = ljnormal(j=jj,y=y,x=x,z=z,nj=nj,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae)
  else if (distr=="st") lv[jj] = ljt(j=jj,nu=nu,y=y,x=x,z=z,nj=nj,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae)
  else if (distr=="ss") lv[jj] = ljs(j=jj,nu=nu,y=y,x=x,z=z,nj=nj,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae)
  else if (distr=="scn") lv[jj] =ljcn(j=jj,nu=nu,y=y,x=x,z=z,nj=nj,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae)
}
  return(lv)
}

################################################################
#Root of a symmetric matrix
################################################################
matrix.sqrt <- function(A)
{
  if (length(A)==1) return(sqrt(A))
  else{
    sva <- svd(A)
    if (min(sva$d)>=0) {
      Asqrt <- sva$u%*%diag(sqrt(sva$d))%*%t(sva$v) # svd e decomposi??o espectral
      if (all(abs(Asqrt%*%Asqrt-A)<1e-4)) return(Asqrt)
      else stop("Matrix square root is not defined/not real")
    }
    else stop("Matrix square root is not defined/not real")
  }
}

#likelihood for envelope
logveroenvel<-function(nj,y,x,z,beta1,sigmae, D1, Gama, lambda, eta, type="MGH"){
  
  m<-length(nj)
  N<-sum(nj)
  p<-dim(x)[2]
  q1<-dim(z)[2]
  
  suma1<-numeric()
  
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
      suma1[j]<-dghyp(y1, object = para, logvalue = TRUE)
    }
  }
  return(suma1)
}
ljnormal <-function(j,y,x,z,nj,beta1,Gammab,Deltab,sigmae){
  c. = -sqrt(2/pi)  
  p= ncol(x);q1=ncol(z)
  y1<-y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))  ]
  x1<-matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
  z1<-matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj) #z1 D1 z1^t + sig2*I_nj ??
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))
}
#
ljt <-function(j,nu,y,x,z,nj,beta1,Gammab,Deltab,sigmae){
  c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  p= ncol(x);q1=ncol(z)
  y1<-y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))  ]
  x1<-matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
  z1<-matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med))
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  #log(2*dmvt(y1,delta = med, sigma = Psi, df = nu,log=F)*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))#veroST1(Psi,Ajj,dj,nu,pp=njj))
  log(2*dtj*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))
}
#
ljs <-function(j,nu,y,x,z,nj,beta1,Gammab,Deltab,sigmae){
  c.=-sqrt(2/pi)*nu/(nu-.5)
  p= ncol(x);q1=ncol(z)
  y1<-y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))  ]
  x1<-matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
  z1<-matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med))
  #f <- function(u) u^(nu - 1)*dmvnorm(y1,med,Psi/u)*pnorm(u^(1/2)*Ajj)
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
  resp <- integrate(Vectorize(f2),0,1)$value
  log(2*nu*resp)
}
#
ljcn <-function(j,nu,y,x,z,nj,beta1,Gammab,Deltab,sigmae){
  c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  p= ncol(x);q1=ncol(z)
  y1<-y[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))  ]
  x1<-matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
  z1<-matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(sqrt(nu[2])*Ajj,0,1)+
           (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
}



##generate data
geradesbalanceado<-function(n,nj,x,z,beta,sigmae,Gama,D1,lambda,eta,type="MGH"){

p<-dim(x)[2]
q1<-dim(z)[2]
 N<-sum(nj)


if (type=="MGH"){theta=c(lambda,eta,eta) }

if (type=="MGHsym"){theta=c(lambda,eta,eta) 
Gama=matrix(0,q1,1) 
}

if (type=="NIG"){theta=c(lambda=-0.5,eta,eta)}

y1=matrix(0,N,1)

for (k in 1:n) {

u1=rgig(n=1,lambda=lambda,chi=eta,psi=eta)
zb=rmvnorm(n=1,mean = rep(0, nrow(D1)),D1)
ze=rmvnorm(n=1,mean = rep(0, nj[k]),sigmae*diag(nj[k]))


b=u1*Gama+t(sqrt(u1)*zb)
e=t(sqrt(u1)*ze)
x1<-matrix(x[(sum(nj[1:k-1])+1) : (sum(nj[1:k])),  ],ncol=p)
z1<-matrix(z[(sum(nj[1:k-1])+1) : (sum(nj[1:k])) ,  ],ncol=q1)
y1[(sum(nj[1:k-1])+1) : (sum(nj[1:k]))  ]=x1%*%beta+z1%*%b+e
}
return(y1)
}


rsmsn.lmm <- function(time1,x1,z1,sigma2,D1,beta,lambda,depStruct="CI",phi=NULL,distr="sn",nu=NULL) {
  if (length(D1)==1 & !is.matrix(D1)) D1=as.matrix(D1)
  q1 = nrow(D1)
  p = length(beta)
  if (ncol(as.matrix(x1))!=p) stop("incompatible dimension of x1/beta")
  if (ncol(as.matrix(z1))!=q1) stop ("incompatible dimension of z1/D1")
  if (length(lambda)!=q1) stop ("incompatible dimension of lambda/D1")
  if (!is.matrix(D1)) stop("D must be a matrix")
  if ((ncol(D1)!=q1)|(nrow(D1)!=q1)) stop ("wrong dimension of D")
  if (length(sigma2)!=1) stop ("wrong dimension of sigma2")
  if (sigma2<=0) stop("sigma2 must be positive")
  Sig <- sigma2*diag((time1))   #alterado para o caso independente: camila
  #
  if (distr=="ssl") distr<-"ss"
  if (!(distr %in% c("sn","st","ss","scn"))) stop("Invalid distribution")
  if (distr=="sn") {ui=1; c.=-sqrt(2/pi)}
  if (distr=="st") {ui=rgamma(1,nu/2,nu/2); c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {ui=rbeta(1,nu,1); c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {ui=ifelse(runif(1)<nu[1],nu[2],1);
                      c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  #if (all(lambda==0)) c.=0
  delta = lambda/as.numeric(sqrt(1+t(lambda)%*%(lambda)))
  Delta = matrix.sqrt(D1)%*%delta
  Gammab = D1 - Delta%*%t(Delta)
  Xi = matrix(x1,ncol=p)
  Zi = matrix(z1,ncol=q1)
  Beta = matrix(beta,ncol=1)
  ti = c.+abs(rnorm(1,0,ui^-.5))
  bi = t(rmvnorm(1,Delta*ti,sigma=ui^(-1)*Gammab))
  Yi = t(rmvnorm(1,Xi%*%Beta+Zi%*%bi,sigma=ui^(-1)*Sig))
  if (all(Xi[,1]==1)) Xi = Xi[,-1]
  if (all(Zi[,1]==1)) Zi = Zi[,-1]
  return(Yi)
}

##########################################################################
## ENVELOPE FUNCTION
##########################################################################
envelope <- function(model,nj,y,x,z,rep=500,level=.99){
  ecdfext <- function(g) ecdf(g)(g)
  vmat <- matrix(nrow=length(nj),ncol=rep)
  empmat <- matrix(nrow=length(nj),ncol=rep)
  
  if(inherits(model,"EM.MGH")){
    
    for(i in 1:rep){
      cat("Monte Carlo",i,"\r")
      ynovo=geradesbalanceado(n=length(nj),nj,x,z,model$beta1,model$sigmae, model$Gama,model$D1,model$lambda, model$eta,type=model$type)
      vmat[,i]=logveroenvel(nj,ynovo,x,z,model$beta1,model$sigmae, model$D1, model$Gama, model$lambda, model$eta, type="MGH")
      empmat[,i]=ecdfext(vmat[,i])
    }
    
    vv <- c(vmat)
    emp <- c(empmat)
    vv <- round(vv,1)
    labels <- sort(unique(vv))
    
    alpha <- (1-level)/2
    minemp <- tapply(emp,vv,function(x) quantile(x,probs=alpha))
    maxemp <- tapply(emp,vv,function(x) quantile(x,probs=(1-alpha)))
    
    
    u=logveroenvel(nj,y,x,z,model$beta1,model$sigmae, model$D1, model$Gama, model$lambda, model$eta, type="MGH")
    empiricau=ecdf(u)
    plot(empiricau,u,xlab="Log density",ylab="ECDF",main="",xlim=range(u))
    text.title <- switch (model$type,
                          "MGH" = "GH - LME",
                          "MGHsym" = "GH sym - LME",
                          "NIG" = "NIG - LME"
    )
    title(text.title)
    polygon(c(labels,labels[length(labels):1]),
            c(minemp,maxemp[length(labels):1]),col=alpha('gray',1.5),border = NA)
    par(new=TRUE)
    plot(empiricau,u,xlab="",ylab="",main="",xlim=range(u))
    return(invisible(NULL))
  }
  if(inherits(model,"SMSN")){
    beta1<-model$est$beta
    sigmae<- model$est$sigma2 
    D1<-model$est$ds
    D1<- matrix(c(D1[1],D1[2],D1[2],D1[3]), 2, 2) 
    D1<-D1%*%D1
    lambda<-model$est$lambda
    if(model$distr =="st") nu <- model$est$nu
    
    N<-sum(nj)
    m<-length(nj)
    p<-dim(x)[2]
    q1<-dim(z)[2]
    
    for(i in 1:rep){
      cat("Monte Carlo",i,"\r")
      ynovo=matrix(0,N,1)
      
      for (j in 1:m){
        x1<-matrix(x[(sum(nj[1:j-1])+1) : (sum(nj[1:j])),  ],ncol=p)
        z1<-matrix(z[(sum(nj[1:j-1])+1) : (sum(nj[1:j])) ,  ],ncol=q1)
        yaux=rsmsn.lmm(time1=nj[j],x1,z1,sigmae,D1,beta1,lambda,depStruct="CI",phi=NULL,distr=model$distr,nu=nu)
        ynovo[(sum(nj[1:j-1])+1) : (sum(nj[1:j]))  ]=yaux
      }
      
      vmat[,i]=logvero.lmm(ynovo,x,z,nj,beta1,sigmae,D1,lambda,distr="sn",nu=null)
      empmat[,i]=ecdfext(vmat[,i])
    } 
    
    vv <- c(vmat)
    emp <- c(empmat)
    vv <- round(vv,1)
    labels <- sort(unique(vv))
    
    alpha <- (1-level)/2
    minemp <- tapply(emp,vv,function(x) quantile(x,probs=alpha))
    maxemp <- tapply(emp,vv,function(x) quantile(x,probs=(1-alpha)))
    
    u=logvero.lmm(y=y,x=x,z=z,ind=nj,beta1,sigmae,D1,lambda,distr=model$distr,nu=nu)
    empiricau=ecdf(u)
    plot(empiricau,u,xlab="Log density",ylab="ECDF",main="",xlim=range(u))
    text.title <- switch (model$distr,
                          "st" = "ST - LME",
                          "sn" = "SN - LME"
    )
    title(text.title)
    polygon(c(labels,labels[length(labels):1]),
            c(minemp,maxemp[length(labels):1]),col=alpha('gray',1.5),border = NA)
    par(new=TRUE)
    plot(empiricau,u,xlab="",ylab="",main="",xlim=range(u))
    return(invisible(NULL))
  }
  stop("model class not recognized.")
}