args=(commandArgs(TRUE))
traits.file			=	args[1] 			# file of intersecting SNPs and their ABFs (OUTtoclump.file from prep_COMET_general_script.R)
R					=	as.numeric(args[2]) # ratio of costs = type II error cost/ type I error cost
pi0					=	as.numeric(args[3]) # Pr(H_0); e.g. 0.99, 0.999
clumped.file		=	args[4] 			# pruned SNPs using ABF measure; plink output that has been gzipped
covariate.file		=	args[5]				# file of covariate values for each SNP
covariate.names		=	unlist(strsplit(args[6],","))				# vector of covariate names
printsig			=	args[7]				# value 1 if significant overlap variants, their marginal statistics, and covariate values are to be output; 0 otherwise
main.output.file	=	args[8]				# output file with COMET results


print(covariate.names)

source("COMET_script.R")

traits <- read.table(gzfile(traits.file),header=TRUE,as.is=TRUE)
CLabf <- read.table(gzfile(clumped.file),header=TRUE,as.is=TRUE)
Cinfo <- read.table(gzfile(covariate.file),header=TRUE,as.is=TRUE)


# match on position, keep clumped snps
keep2 <- match(CLabf$chrBP,Cinfo$chrBP)
Ckeep2 <- Cinfo[keep2[!is.na(keep2)],]
keepT <- match(Ckeep2$chrBP,traits$chrBP)
tabf <- traits[keepT[!is.na(keepT)],]


# check that unpruned and clumped sets of proportions are similar
keepC <- match(traits$chrBP,Cinfo$chrBP)
Ctraits <- Cinfo[keepC[!is.na(keepC)],]

print(apply(Ctraits[,-(1:2)],2,mean)) # unpruned proportions
print(apply(Ckeep2[,-(1:2)],2,mean)) # uclumped proportions


out <- fitmodels.fn(t1_abf=tabf$t1_abf,t2_abf=tabf$t2_abf,R=R,pi0=pi0,Cinfo=Ckeep2,covnames=covariate.names,printsig=printsig)
#print(out)

write.table(out$models$out1,main.output.file,quote=FALSE,sep="\t",append=FALSE)
write.table("\n",main.output.file,quote=FALSE,sep="\t",append=TRUE)
write.table(out$models$out2,main.output.file,quote=FALSE,sep="\t",append=TRUE)
write.table("\n",main.output.file,quote=FALSE,sep="\t",append=TRUE)
write.table(out$models$out12,main.output.file,quote=FALSE,sep="\t",append=TRUE)



if(printsig==1) write.table(out$overlap,paste(main.output.file,"-overlap-details.txt",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")


q(save="no");
