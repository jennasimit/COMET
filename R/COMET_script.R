### This script is called by COMET_general_script.R

#### define functions


##
inv.logit <- function(x) exp(x)/(1+exp(x))

##
logit <- function(p) log(p/(1-p))

##
# 2-step inter-rater models
fitmodels.fn <- function(t1_abf,t2_abf,R=R,pi0=pi0,Cinfo,covnames,printsig)	{ 

PO <- pi0/(1-pi0)
theta <- PO/R

sigabf1 <- 1*(t1_abf >= theta)
sigabf2 <- 1*(t2_abf >= theta)

 # e.g. covnames <- c("Q1","Q2","Q3","Q4","Q5","Q6")
 Ckeep <- Cinfo[,c("chrBP",covnames)]
 
 # marginal models
 mod1 <- glm(sigabf1~.-chrBP,Ckeep,family=binomial(link = "logit"))
 mod2 <- glm(sigabf2~.-chrBP,Ckeep,family=binomial(link = "logit"))

# calculate offset term, eta.assoc
c1 <- coefficients(mod1) 
c2 <- coefficients(mod2)

 a1 <- c1[1]
 a2 <- c2[1]

K <- dim(Ckeep)[2]
for(j in 2:K) {
 a1 <- a1 + c1[j]*Ckeep[,j]
 a2 <- a2 + c2[j]*Ckeep[,j]
  }

 
 p1 <- inv.logit(a1)
 p2 <- inv.logit(a2)
  
  
 sigabf12 <- sigabf1*sigabf2
 eta.assoc <- logit(p1*p2)

# before fitting overlap model, remove any covariate categories that contain 0 overlap variants 
check <- apply(Ckeep[which(sigabf12==1),-1],2,mean)
rmcov <- covnames[which(check==0)]
Ckeep2 <- Ckeep
if(length(rmcov) > 0 & length(rmcov)< K) { keepcov <- setdiff(covnames,rmcov);  Ckeep2 <- Ckeep[,c("SNP",keepcov)] }

if(length(rmcov)==(K-1))	{
 output <- summary(mod12.full)$coefficients
 output[covnames,1:3] <- 0
 output[covnames,4] <- 1
  					}

if(length(rmcov)<K)	{  
mod12.assoc <- glm(sigabf12~.-chrBP,Ckeep2,offset=eta.assoc,family=binomial(link = "logit"))

output <- summary(mod12.assoc)$coefficients
 
# fit a reduced overlap model, removing any covariates that are poorly estimated (i.e. have stderr > 3)
 Qrm <- which(summary(mod12.assoc)$coefficients[,2] > 3)
 K2 <- dim(Ckeep2)[2]
 if(length(Qrm)>0 & length(Qrm)< K2)	{
 ccred <- setdiff(row.names(summary(mod12.assoc)$coefficients),c(names(Qrm),"(Intercept)"))
 Qinfored <- Ckeep2[,ccred]
 mod12.full <- mod12.assoc 
 mod12.assoc <- glm(sigabf12~.,data=Qinfored,offset=eta.assoc,family=binomial(link = "logit"))
 
 output <- summary(mod12.full)$coefficients
 output[row.names(summary(mod12.assoc)$coefficients),] <- summary(mod12.assoc)$coefficients
 output[names(Qrm),1:3] <- 0
 output[names(Qrm),4] <- 1
  					}
 
 if(length(Qrm)==K2)	{
 output <- summary(mod12.assoc)$coefficients
 output[-1,1:3] <- 0
 output[-1,4] <- 1
  					}
 
 }
  
 out1 <- summary(mod1)
 out2 <- summary(mod2)
 out12 <- output
 
 Px <- function(x) 1-pnorm(x)
 
 # determine the number of associated (marginal: c1, c2; overlap: c12) variants in each covariate category
 c1 <- c(NA,apply(Ckeep[sigabf1==1,-1],2,sum)) # NA is given for the "number of associated variants" for the intercept
 c2 <- c(NA,apply(Ckeep[sigabf2==1,-1],2,sum))
 c12 <- c(NA,apply(Ckeep2[sigabf12==1,-1],2,sum))
 
 pv1 <- c(coef(out1)[1,4], Px(coef(out1)[-1,3]))
 pv2 <- c(coef(out2)[1,4], Px(coef(out2)[-1,3]))
 pv12 <- c(out12[1,4], Px(out12[-1,3]))


 
 o1 <- data.frame(pv1,coef(out1)[,1:2],c1)
 names(o1) <- c("P","Effect","StdErr","Count")
 o2 <- data.frame(pv2,coef(out2)[,1:2],c2)
 names(o2) <- c("P","Effect","StdErr","Count")
 o12 <- data.frame(pv12,out12[,1:2],c12)
 names(o12) <- c("P","Effect","StdErr","Count")


  models <-  list(out1=t(o1),out2=t(o2),out12=t(o12))
 
 if(printsig==0) return(models=models)
 				
 
 if(printsig==1) {
  overlap <- data.frame(Cinfo[which(sigabf12==1),c("SNP","chrBP",covariate.names)],t1_L10abf=log10(t1_abf[which(sigabf12==1)]),t2_L10abf=log10(t2_abf[which(sigabf12==1)]))
  return(list(models=models,overlap=overlap))
 				}
 				
 				
       }

