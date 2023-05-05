# Full conditional distribution for beta
atualizarBETA <- function(mu,sigma,y,x,u,tau2){
  s.inv <- chol2inv(chol(sigma))
  mcov  <- chol2inv(chol(s.inv+(1/tau2)*(t(x)%*%x)))
  v     <- y-u
  media <- ((1/tau2)*t(v)%*%x+mu%*%s.inv)%*%mcov
  beta  <- rmvnorm(1,media,mcov)
  return(beta)
}
# Full conditional distribution for tau2
atualizarTAU2<-function(alpha0,beta0,n,l,y,x,beta,u){
  beta  <- beta0 + sum((y-x%*%beta-u)^2)/2
  alpha <- alpha0 + n*l/2
  tau2  <- 1/rgamma(1,alpha,beta)
  return(tau2)
}
# Full conditional distribution for inefficiencies u
condicionalU<-function(u,eta,sigma2,y,x,beta,tau2){
  priori <- -log(u)-(log(u)-eta)^2/(2*sigma2)
  verossi<- -(y-x%*%beta-u)^2/(2*tau2)
  funcao <- priori+verossi
  return(funcao)
}
# Metropolis-Hasting for inefficiencies u
atualizarU<-function(u,eta,sigma2,y,x,beta,tau2){
  valoratual   <- u
  valorproposto<- 0.01+exp(rnorm(1,log(valoratual-0.01),1))
  candidato    <- exp(condicionalU(valorproposto,eta,sigma2,y,x,beta,tau2)-condicionalU(valoratual,eta,sigma2,y,x,beta,tau2)+log(valorproposto-0.01)-log(valoratual-0.01))
  chanceaceitar <- min(1,candidato)
  Ufinal        <- ifelse(runif(1)<chanceaceitar,valorproposto,valoratual)  
  return(Ufinal)
}
# Full conditional distribution for sigma2
atualizarSIGMA2<-function(alpha1,beta1,n,l,u,eta){
  w      <- log(u)
  alpha  <- alpha1 + 0.5*(n*l)
  beta   <- beta1 + sum((w-eta)^2)/2
  sigma2 <- 1/rgamma(1,alpha,beta)
  return(sigma2)
}
# Full conditional distribution for omega
atualizarW2<-function(alpha2,beta2,n,l,quad.sum){
  alpha <- alpha2 + n*l
  beta  <- beta2  + quad.sum
  w2    <- riwish(alpha, beta)
  return(w2)
}
# Full conditional for the latent variable in the degrees of freedom
atualizarV<-function(nu,quad.nmult,l){
  alpha <- c(rep((nu+1)/2,l))
  beta  <- nu/2 + quad.nmult
  v     <- rgamma(l,alpha,beta)
  return(v)
}
# Full conditional for the degrees of freedom 
condicionalNU1<-function(nu,v,n,l){
  d       <- 4/(1+sqrt(2))
  priori  <- log(nu)-3*log(nu+d)
  verossi <- (n*l*nu/2)*log(nu/2)-n*l*log(gamma(nu/2))+(nu/2-1)*sum(log(v))-nu/2*sum(v)
  funcao  <- priori+verossi
  return(funcao)
}
# Metropolis-Hasting for the degrees of freedom 
atualizarNU1<-function(nu,v,n,l,clap,clap.aux,M0,t){
  valoratual<-nu
  valorproposto<-rtnorm(1,valoratual,sqrt(clap*clap.aux),lower=2,upper=40)
  candidato<-exp(condicionalNU1(valorproposto,v,n,l)-condicionalNU1(valoratual,v,n,l)-dtnorm(valorproposto,valoratual,sqrt(clap*clap.aux),lower=2,upper=40,log=TRUE)+dtnorm(valoratual,valorproposto,sqrt(clap*clap.aux),lower=2,upper=40,log=TRUE))
  chanceaceitar<-min(1,candidato)
  NUfinal <- ifelse(runif(1)<chanceaceitar,valorproposto,valoratual)
  gama1<- 1/t^0.5
  gama2<- 1*gama1
  termometro<- exp(log(clap)+gama2*(chanceaceitar-0.234))
  termometro.aux<- clap.aux+gama1*((NUfinal-M0)^2-clap.aux)
  p.auxiliar<-M0+gama2*(NUfinal-M0)
  return(c(NUfinal,termometro,termometro.aux,p.auxiliar))
}