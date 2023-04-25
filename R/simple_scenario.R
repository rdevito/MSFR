rm(list=ls())
library(matlab)
library(R.matlab)
library(MASS)
library(GPArotation)
library(psych)
#install.packages('R.matlab', dependencies=TRUE, repos='http://cran.rstudio.com/')

#source('/Users/roberta/Dropbox/CrossStudyBFR/Rcode/MSFA_EM_beta/MSFA_R_beta.R')
#source("~/Dropbox (HMS)/Harvard/EMCrossStudyBFR/Rcode/MSFA_EM_beta/MSFR.R")
source('/Users/roberta/Dropbox/CrossStudyBFR/Rcode/MSFA_EM_beta/MSFR.R')

#############################
set.seed(202110)
###generate the data: really simple
S             <- 2
p             <- 20
p_b           <- 2
k             <- 3
j_s           <- rep(1, S)
n_s           <- rep(500, S)
#n_s          <- c(10, 15, 12, 14)
#theta        <- rep(0, length = p)
PH            <- as.vector(zeros(p, k))
noZEROc       <- (p / 3) * k
studyc        <- runif(noZEROc, 0.6, 1)
sign          <- sample(x = length(studyc), size = (length(studyc) / 2))
studyc[sign]  <- studyc[sign] * (-1)
positionc     <- sample(x = k * p, size = length(studyc))
PH[positionc] <- studyc
Phi           <- matrix(PH, p, k)
Beta          <- matrix(c(rep(-1,p/2),rep(1,p/2),
                          rep(0,p/4),rep(1,3*p/4)), p, p_b)
L             <- noZERO <- study <- position <- Lambda_s <- 
  Psi_s <- Sigma_s <- X_s <- I_j <-B_s <- Mean_X_s<-list()

for(s in 1:S){
  L[[s]]                <- as.vector(zeros(p, j_s[s]))
  noZERO[[s]]           <- (p / 3) * j_s[s]
  study[[s]]            <- runif(noZERO[[s]], -1, 1)
  position[[s]]         <- sample(x = p * j_s[s], size = length(study[[s]]))
  L[[s]][position[[s]]] <- study[[s]]
  Lambda_s[[s]]         <- matrix(L[[s]], p, j_s[s])
  Psi_s[[s]]            <- diag(runif(p, 0, 1), p)
  Sigma_s[[s]]          <- tcrossprod(Phi) + tcrossprod(Lambda_s[[s]]) + Psi_s[[s]]
  B_s[[s]]              <- matrix(runif(n_s[s]*p_b,0, 1),n_s[s],p_b)
  Mean_X_s[[s]]         <- B_s[[s]]%*%t(Beta)
  X_s[[s]]              <- do.call("rbind",plyr::llply(.data = 1:n_s[[s]], .fun = function(x){
    mvrnorm(1, Mean_X_s[[s]][x,], Sigma_s[[s]])
  }, .parallel = FALSE))
}

test <- start_msfa(X_s, B_s, p_b, k, j_s, constraint = "block_lower2", method = "adhoc")
test$Lambda_s[[1]]
test$beta

####Checking initialization of betas ####
B <- Reduce(rbind, B_s)
X <- Reduce(rbind, X_s)
B <- Reduce(rbind, B_s)
X <- Reduce(rbind, X_s)

#With R function
#this function works because X and B don't have names and 
#we don't consider categorical data
fm1 <- lm(X ~ 0+B) 
beta_fm1 = t(fm1$coefficients)

#Manual
beta_mle = t(solve(t(B) %*% (B)) %*% t(B) %*% X)
round(beta_mle-beta_fm1,10) #all 0's :)

t1<-Sys.time()
EM_beta <- ecm_msfa(X_s, B_s, start=test, trace = FALSE)
t2<-Sys.time()
t2-t1

RV=data.frame(Beta=FactoMineR::coeffRV(EM_beta$beta,Beta)$rv,
              Phi=FactoMineR::coeffRV(EM_beta$Phi,Phi)$rv,
              Lambda1=FactoMineR::coeffRV(EM_beta$Lambda_s[[1]],Lambda_s[[1]])$rv,
              Lambda2=FactoMineR::coeffRV(EM_beta$Lambda_s[[2]],Lambda_s[[2]])$rv,
              Psi1=FactoMineR::coeffRV(diag(EM_beta$psi_s[[1]]),Psi_s[[1]])$rv,
              Psi2=FactoMineR::coeffRV(diag(EM_beta$psi_s[[2]]),Psi_s[[2]])$rv,
              Sigma1=FactoMineR::coeffRV(tcrossprod(EM_beta$Phi) + 
                                           tcrossprod(EM_beta$Lambda_s[[1]]) + 
                                           diag(EM_beta$psi_s[[1]]),
                                         tcrossprod(Phi) + tcrossprod(Lambda_s[[1]]) + Psi_s[[1]])$rv,
              Sigma2=FactoMineR::coeffRV(tcrossprod(EM_beta$Phi) + 
                                           tcrossprod(EM_beta$Lambda_s[[2]]) + 
                                           diag(EM_beta$psi_s[[2]]),
                                         tcrossprod(Phi) + tcrossprod(Lambda_s[[2]]) + Psi_s[[2]])$rv)
#0.9900578 0.9963294 0.9667498 0.9843553 0.5390422 0.6387878 0.9937982 0.9957357 #MSFA
