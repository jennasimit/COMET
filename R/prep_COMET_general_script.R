args=(commandArgs(TRUE))
trait1.file		=	args[1] 
trait2.file		=	args[2] 
R			=	as.numeric(args[3])
pi0			=	as.numeric(args[4])
abf.calc	=	args[5] # TRUE or FALSE
OUTtoclump.file	=	args[6]
cc				=	args[7] # TRUE (case-control) or FALSE (QT)
trait1SE			=	as.numeric(args[8])		# set to 1 for all scenarios, unless doing QT analysis and trait is not standard normal nor inverse-rank normalised trait, in that case this use an estimate of the trait's standard error 
trait2SE			=	as.numeric(args[9])		# set to 1 for all scenarios, unless doing QT analysis and trait is not standard normal nor inverse-rank normalised trait, in that case this use an estimate of the trait's standard error 
chr				=	as.numeric(args[10])


# The output file OUTtoclump.file is input into clump-using-plink.sh


# run for each chromosome

source("prep_COMET.R")

t1 <- read.table(trait1.file,header=TRUE)
t2 <- read.table(trait2.file,header=TRUE)


# if already have ABFs for each trait (column names: t1_abf, t2_abf), then set abf.calc=FALSE and this function will output the data frame of intersecting SNPs 
# if have summary statistics, then set abf.calc=TRUE and this function will output the data frame of intersecting SNPs and calculate the marginal ABFs 
# first four columns: SNP	BP	REF	ALT

abf.data <- data.abf.fn(t1=t1,t2=t2,cc=cc,abf.calc=abf.calc,traitSE1=trait1SE,traitSE2=trait2SE)

abf.measure <- clump.measure.fn(traits=abf.data,R=R,pi0=pi0,chr=chr)
write.table(abf.measure,file=OUTtoclump.file,row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")


q(save="no");
