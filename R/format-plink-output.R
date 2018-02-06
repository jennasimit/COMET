args=(commandArgs(TRUE))
clumped.file		=	args[1] 			# pruned SNPs using joint measure; plink output 
chr					=	as.numeric(args[2])


chrbp.fn <- function(x,chr,j) {
 if(is.na(x$CHR[j])) x$BP[j] <- x$SNP[j]
 else x$BP[j] <- paste(chr,":",x$BP[j],sep="")
 return(x$BP[j])
 								}

CLtraits <- read.table(clumped.file,header=TRUE,as.is=TRUE)
m <- dim(CLtraits)[1]
m.mat <- matrix(1:m,nrow=m)
chrbp <- apply(m.mat,1,chrbp.fn,x=CLtraits,chr=chr)
CLtraits$BP <- chrbp
names(CLtraits)[4] <- "chrBP"	

write.table(CLtraits,clumped.file,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
