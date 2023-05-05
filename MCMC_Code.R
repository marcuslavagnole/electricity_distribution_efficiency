# Load libraries
library(mvtnorm)
library(MCMCpack)
library(coda)
library(msm)
#library(foreign)
#library(ggplot2)
#library(reshape2)
#library(openxlsx)
#library(plyr)

# Set working directory to file location and load utils
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Full_Conditionals.R")
load("base.RData")

# Set appropriately the dependent and indepent variables 
CC = NULL
for(i in 1:(length(base$Ano)/14)){
  CC = cbind(CC,base$Custo[((i-1)*14+1):(i*14)])
}
CC = log(CC)
QQ = NULL
for(i in 1:(length(base$Ano)/14)){
  QQ = cbind(QQ,base$UC[((i-1)*14+1):(i*14)])
}
QQ = log(QQ)
PF = NULL
for(i in 1:(length(base$Ano)/14)){
  PF = cbind(PF,base$Mercado[((i-1)*14+1):(i*14)])
}
PF = log(PF)
LF = NULL
for(i in 1:(length(base$Ano)/14)){
  LF = cbind(LF,base$Rede[((i-1)*14+1):(i*14)])
}
LF = log(LF)

# Define panel dimension
l<- dim(CC)[1]
n<- dim(CC)[2]

# Define number of iterations
NN<-100000

# Create auxiliary objects
theta2 <- matrix(NA,(2*n),(l+1))
omega2 <- array(dim=c(2,2,NN))

theta1 <- array(dim=c(l,n,NN))
sigma2 <- matrix(NA,NN,1)

beta   <- matrix(NA,NN,4)
tau2   <- matrix(NA,NN,1)

a<- list() ; R<- list()
m<- list() ; C<- list()

nu <- matrix(NA,NN,1)
nu1 <- matrix(NA,NN,1)

# Set the initial values
m[[1]]<- c(0,0)
C[[1]]<- diag(c(1,1))

theta1[,,1] <- matrix(rep(rep(0.1,n),l),ncol=n,byrow=TRUE)
omega2[,,1] <- c(1,0,0,1)
sigma2[1,1] <- 1
tau2[1,1]   <- 1
beta[1,]    <- c(-1,0.6,0.2,0.15)

v <- matrix(1,14,60)
nu1[1,1] <- 4

# Create auxiliary objects for the adaptative MH
rmw_nu1     <- matrix(NA,NN,3)
rmw_nu1[1,] <- c(0.8,1,0)

# MCMC
pb <- txtProgressBar(min = 2, max = NN, style = 2)
for (i in 2:NN){
  omega.aux <- matrix(0,2,2)
  v.aux     <- matrix(NA,l,n)
  #Forward Filtering
  V<- sigma2[i-1,1]
  for (s in 1:n){
    for (k in 2:(l+1)){
      W<- omega2[,,i-1]*(1/v[(k-1),s])
      # Prior
      G <- matrix(c(1,1,0,1),2,byrow=TRUE) #2x2
      a[[k]] <- G%*%m[[k-1]]          #2x1      
      R[[k]] <- G%*%C[[k-1]]%*%t(G)+W #2x2
      # One-step forecast
      F1  <- matrix(c(1,0),2,1)     #2x1
      f   <- t(F1)%*%a[[k]]         #1x1
      Q   <- t(F1)%*%R[[k]]%*%F1+V  #1x1
      #Posterior
      A      <- R[[k]]%*%F1%*%solve(Q)                   #2x1
      m[[k]] <- a[[k]]+A%*%(log(theta1[k-1,s,i-1])-c(f)) #2x1
      C[[k]] <- R[[k]]-A%*%Q%*%t(A)                      #2x2
    }
    #Backward Sampling
    C[[(l+1)]][lower.tri(C[[(l+1)]])] <- C[[(l+1)]][upper.tri(C[[(l+1)]])]
    theta2[(2*(s-1)+1):(2*s),(l+1)]<- rmvnorm(1,m[[(l+1)]],C[[(l+1)]])
    for (j in 1:l){
      B <- C[[(l+1-j)]]%*%t(G)%*%solve(R[[(l+2-j)]])
      h <- m[[(l+1-j)]]+B%*%(theta2[(2*(s-1)+1):(2*s),(l+2-j)]-a[[(l+2-j)]])
      H <- C[[(l+1-j)]]-B%*%R[[(l+2-j)]]%*%t(B)
      H[lower.tri(H)] <- H[upper.tri(H)]
      theta2[(2*(s-1)+1):(2*s),(l+1-j)]<- rmvnorm(1,h,H)
      
      theta1[(l+1-j),s,i] <- atualizarU(theta1[(l+1-j),s,i-1],(t(c(1,0))%*%theta2[(2*(s-1)+1):(2*s),(l+2-j)]),sigma2[i-1,1],CC[(l+1-j),s],c(1,QQ[(l+1-j),s],PF[(l+1-j),s],LF[(l+1-j),s]),beta[i-1,],tau2[i-1,1])

      v.aux[(l+1-j),s]    <- 0.5*t(theta2[(2*(s-1)+1):(2*s),(l+2-j)]-G%*%theta2[(2*(s-1)+1):(2*s),(l+1-j)])%*%chol2inv(chol(omega2[,,i-1]))%*%(theta2[(2*(s-1)+1):(2*s),(l+2-j)]-G%*%theta2[(2*(s-1)+1):(2*s),(l+1-j)])
    }
    v[,s]    <- atualizarV(nu1[i-1,1],v.aux[,s],l)
    omega.aux <- omega.aux + (rbind(v[,s],v[,s])*(theta2[(2*(s-1)+1):(2*s),-1]-G%*%theta2[(2*(s-1)+1):(2*s),-(l+1)]))%*%t(theta2[(2*(s-1)+1):(2*s),-1]-G%*%theta2[(2*(s-1)+1):(2*s),-(l+1)])
  }
  #Parameters
  sigma2[i,1]  <- atualizarSIGMA2(0.1,0.1,n,l,c(theta1[,,i]),c(t(theta2[seq(1,(2*n),2),-1]))) 
  
  omega2[,,i] <- atualizarW2(3,diag(c(1,1)),n,l,omega.aux)
  nu1.aux     <- atualizarNU1(nu1[i-1,1],c(v),n,l,rmw_nu1[i-1,1],rmw_nu1[i-1,2],rmw_nu1[i-1,3],i)
  nu1[i,1]    <- nu1.aux[1] ; rmw_nu1[i,] <- nu1.aux[2:4]
  
  beta[i,]    <- atualizarBETA(c(0,0,0,0),diag(10,4),c(CC),cbind(rep(1,n*l),c(QQ),c(PF),c(LF)),c(theta1[,,i]),tau2[i-1,1])
  tau2[i,1]   <- atualizarTAU2(0.1,0.1,n,l,c(CC),cbind(rep(1,n*l),c(QQ),c(PF),c(LF)),beta[i,],c(theta1[,,i]))
  gc()
  setTxtProgressBar(pb, i)
} 
