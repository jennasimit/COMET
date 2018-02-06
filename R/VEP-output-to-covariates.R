args=(commandArgs(TRUE))	
fname	=	as.character(args[1]) # vep.output
fout	=	as.character(args[2]) # covariate output

# Here, VEP was run with the settings: --terms so --check_existing --regulatory 

vep <- read.table(gzfile(fname),as.is=TRUE,fill=TRUE,col.names=c("Uploaded_variation","Location","Allele","Gene","Feature","Feature_type","Consequence","cDNA_position","CDS_position","Protein_position","Amino_acids","Codons","Existing_variation","Extra")) 


#Q1:SNP in a transcribed but not translated region
#Q2:SNP in a translated region but does not change the amino acid
#Q3:SNP is potentially deleterious (might disrupt protein)
#Q4:SNP is potentially loss of function (might cause loss of function of protein) (defined as in lofthee  of MacArthur et al. 2012)
#Q5:SNP is in a potentially regulatory or regulatory region
#Q6:SNP is intergenic

q1 <- c("mature_miRNA_variant", "non_coding_transcript_exon_variant", "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant")
q2 <- c("stop_retained_variant","synonymous_variant","incomplete_terminal_codon_variant")
q3 <- c("inframe_insertion", "inframe_deletion", "missense_variant","initiator_codon_variant","splice_region_variant","stop_lost")
q4 <- c("stop_gained","frameshift_variant","splice_donor_variant","splice_acceptor_variant","start_lost")
q5 <- c("5_prime_UTR_variant","3_prime_UTR_variant","TF_binding_site_variant","regulatory_region_ablation","regulatory_region_amplification","regulatory_region_variant")
q6 <- c("intergenic_variant","upstream_gene_variant","downstream_gene_variant")
q7 <- c("splice_donor_variant","splice_acceptor_variant","stop_gained","frameshift_variant") # lof1
q8 <- c("splice_donor_variant","splice_acceptor_variant","stop_gained","frameshift_variant","stop_lost","initiator_codon_variant") # lof2
q9 <- c("stop_lost","initiator_codon_variant") # lof2-lof1
q10 <- "missense_variant"

question.fn <- function(vepout,k)       {
 out <- numeric(10)
 if( sum(charmatch(unlist(strsplit(as.character(vepout[k]),",")),q1,nomatch=0))>0)  out[1] <- 1
 if( sum(charmatch(unlist(strsplit(as.character(vepout[k]),",")),q4,nomatch=0))>0 & out[1]==0)  out[4] <- 1
 if( sum(charmatch(unlist(strsplit(as.character(vepout[k]),",")),q3,nomatch=0))>0 & out[1]==0 & out[4]==0)  out[3] <- 1
 if( sum(charmatch(unlist(strsplit(as.character(vepout[k]),",")),q2,nomatch=0))>0 & out[1]==0 & out[3]==0 & out[4]==0)  out[2] <- 1
 if( sum(charmatch(unlist(strsplit(as.character(vepout[k]),",")),q5,nomatch=0))>0 & out[1]==0 & out[3]==0 & out[4]==0)  out[5] <- 1
 if( sum(charmatch(unlist(strsplit(as.character(vepout[k]),",")),q6,nomatch=0))>0 & out[1]==0 & out[3]==0 & out[4]==0  & out[2]==0)  out[6] <- 1

 if( sum(charmatch(unlist(strsplit(as.character(vepout[k]),",")),q7,nomatch=0))>0 & out[1]==0)  out[7] <- 1

 if( sum(charmatch(unlist(strsplit(as.character(vepout[k]),",")),q8,nomatch=0))>0 & out[1]==0)  out[8] <- 1
    if( sum(charmatch(unlist(strsplit(as.character(vepout[k]),",")),q9,nomatch=0))>0 & out[1]==0)  out[9] <- 1
 if( sum(charmatch(unlist(strsplit(as.character(vepout[k]),",")),q10,nomatch=0))>0 & out[1]==0)  out[10] <- 1
 return(out)
                                                                        }


combine.fn <- function(ordvep) 1*(sum(ordvep)>0)


keeprs.fn <- function(x,j) 	{
 check <- unlist(strsplit(x[j,"Existing_variation"],","))
 rsid <- grep("rs",check,value=TRUE)
 m <- length(rsid)
 if(m==0) {out <- x[j,]; out[,"Existing_variation"]=NA}
 if(m==1) {out <- x[j,]; out[,"Existing_variation"]=rsid}
 if(m>1) {
  out <- NULL
  for(k in 1:m) {
   out <- rbind(out,x[j,])
   out[k,"Existing_variation"] <- rsid[k]
  		  		}
  		  }		
 return(out)
							}

trs <- vep[grep("rs",vep[,"Uploaded_variation"]),]  
trs[,"Existing_variation"] <- trs[,"Uploaded_variation"] # rsid in dataset  # set 1

trs.more <- vep[grep("rs",vep[,"Uploaded_variation"],invert=TRUE),]  # no rsid in dataset
trs2 <- trs.more[grep(",",trs.more[,"Existing_variation"],invert=TRUE),] # no rsid in dataset and one co-location given
trs2a <- trs2[grep("rs",trs2[,"Existing_variation"]),] # no rsid in dataset and one co-location given - only rsid co-locations here  # set 2

trs.more2 <- trs.more[grep(",",trs.more[,"Existing_variation"]),] # no rsid in dataset and multiple co-locations given

m <- dim(trs.more2)[1]
trs.more3a <- NULL

if(m>0) {
m.mat <- matrix(1:m,nrow=m)  
trs.more3 <- apply(m.mat,1,keeprs.fn,x=trs.more2) 
trs.more3a <- do.call("rbind",trs.more3)  # set 3, no rsid in dataset and multiple co-locations given, but keep only rsid or NA if no rsid co-location
}

vepF <- rbind(trs,trs2a,trs.more3a)	

# covariate assignment 		
					
m <- length(unique(vepF[,"Existing_variation"]))
tmp <- unique(vepF[,"Existing_variation"])
uni <- as.character(levels(tmp))[tmp]# length m
vepsnp <- vepF[,"Existing_variation"] # snp ids
veppos <- vepF[,"Location"]
uniF <- numeric(m)# length m
vepsnpF <- vepF[,"Consequence"] # consequences 

n <- dim(vepF)[1]
k.mat <- matrix(1:n,nrow=n)
Q.out <- apply(k.mat,1,question.fn,vepout=vepsnpF)
Q.out <- t(Q.out) # columns are Q1-Q6


 Lvec <- 1:length(vepsnp)
 vepall <- data.frame(vepsnp=vepsnp,veppos=veppos,Q.out=Q.out)
 ordsnp <- order(vepall[,"veppos"])
 ordvepall <- vepall[ordsnp,]
 same <- match(ordvepall[,"veppos"],ordvepall[,"veppos"])
 usame <- unique(same) # indices of unique snps in same order as ordvepall

 F2<-c()
 K <- 2+dim(Q.out)[2]
 for(s in usame) { k<- Lvec[same==s]; F2 <- rbind(F2,apply(ordvepall[k,3:K],2,combine.fn))}

out <- data.frame(snp=ordvepall[usame,1:2],F=F2)

names(out) <- c("SNP","chr.BP",paste("Q",1:12,sep="")

write.table(out,fout,row.names=FALSE,col.names=TRUE,quote=FALSE)
					

					
