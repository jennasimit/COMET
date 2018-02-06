
# abf calculation

# case-control data
 abfcc.fn <- function(i,beta,se){
   W <- .04
   V <- se[i]^2
   z <- beta[i]/sqrt(V)
   t1 <- V/(V+W)
   t2 <- W/(V+W)
   abf <- sqrt(t1)*exp(t2*(z^2)/2)
   return(abf)
  }

# QT data
# if standard normal or inverse-rank normalised trait, set traitSE=1. Otherwise, need to provide traitSE=standard error of trait
  abfqt.fn <- function(i,beta,se,traitSE=1){
   V <- se[i]^2
   z <- beta[i]/sqrt(V)
   
  W <- .02*traitSE^2
      
   t1 <- V/(V+W)
   t2 <- W/(V+W)
   abf <- sqrt(t1)*exp(t2*(z^2)/2)
   return(abf)
  }



######

data.abf.fn <- function(t1,t2,cc,abf.calc=TRUE,traitSE1=1,traitSE2=1)	{
# if already have ABFs for each trait (column name: "ABF" for each trait), then set abf.calc=FALSE and this function will output the data frame of intersecting SNPs 
# if have summary statistics, then set abf.calc=TRUE and this function will output the data frame of intersecting SNPs and calculate the marginal ABFs 
# first four columns: SNP	BP	REF	ALT

# only consider SNPs present in both data sets by chr position
intsnp <- intersect(t1$BP,t2$BP)
trait1 <- t1[match(intsnp,t1$BP),]
trait2 <- t2[match(intsnp,t2$BP),]
m <- dim(trait1)[1]

if(abf.calc=="yes") {

if(cc=="yes")	{
attach(trait1)
nmat <- matrix(1:m,nrow=m)
t1_abf <- apply(nmat,1,abfcc.fn,beta=Effect,se=StdErr) # 2 by m output
detach(trait1)							
rm(t1)

attach(trait2)
nmat <- matrix(1:m,nrow=m)
t2_abf <- apply(nmat,1,abfcc.fn,beta=Effect,se=StdErr) # 2 by m output
detach(trait2)				
rm(t2)

		}
		
if(cc=="no")	{
attach(trait1)
nmat <- matrix(1:m,nrow=m)
t1_abf <- apply(nmat,1,abfqt.fn,beta=Effect,se=StdErr,traitSE=traitSE1) # 2 by m output
detach(trait1)							
rm(t1)

attach(trait2)
nmat <- matrix(1:m,nrow=m)
t2_abf <- apply(nmat,1,abfqt.fn,beta=Effect,se=StdErr,traitSE=traitSE2) # 2 by m output
detach(trait2)				
rm(t2)
		}


}

if(abf.calc=="no") {
t1_abf <- trait1$ABF
t2_abf <- trait2$ABF
					}
					
traits <- data.frame(trait1[,c("SNP","BP","REF","ALT")],t1_abf=t1_abf,t2_abf=t2_abf)

return(traits)

								
								}
							
#######

clump.measure.fn <- function(traits,R,pi0,chr)	{
# traits data/frame has columns: SNP	BP	REF	ALT ... t1_abf t2_abf

attach(traits)
								
PO <- pi0/(1-pi0) # prior odds of no association
theta <- PO/R

M <- max(t1_abf,t2_abf,na.rm=TRUE)
sig1 <- 1*(t1_abf > theta)
sig2 <- 1*(t2_abf > theta)


sigabfSome <- sig1+sig2
tabf <- pmax(t1_abf,t2_abf) + sigabfSome*M  # use most sig value of all traits at the SNP and if more than one are sig then boost value to increase chance of not pruning out;
# this accounts for the scenarios where both traits are marginally sig at one SNP, and at a SNP in LD with it, only one trait meets sig and has a very large sig value -->
# taking the most sig value at each SNP would mean we keep the SNP with only one trait assoc and prune out the one that is sig for both traits
traits$abf <- tabf
traits$neg.abf <- -tabf
		
detach(traits)
traits$chrBP <- paste(chr,":",traits$BP,sep="")

return(traits)

	}							



