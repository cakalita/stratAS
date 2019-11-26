library('VGAM')
library(plyr)
library("optparse")
library(data.table)

option_list = list(
	make_option("--inp_counts", action="store", default=NA, type='character',
              help="Path to file containing allelic counts for this individual [required]"),
	make_option("--inp_cnv", action="store", default=NA, type='character',
              help="Path to file containing CNV boundaries for this individual [optional]"),	
	make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
	make_option("--min_cov", action="store", default=5 , type='integer',
              help="Minimum number of REF reads and ALT reads to include this site. [default: %default]"),
    make_option("--multi_ind", action="store_true", default=TRUE ,
              help="Does the counts file contain counts from multiple individuals. [default: %default]")              
)
opt = parse_args(OptionParser(option_list=option_list))
DO.multi_ind = opt$multi_ind

vcf = read.table( opt$inp_counts , as.is=T , header = FALSE)
#columns 1:3 chr,pos,RSID,ref then every 4 are sample,hap,ref,alt
if (DO.multi_ind) {
vcf_list <- lapply(seq(5, ncol(vcf), by=4), function(i) {
	#i<-seq(5, ncol(vcf), by=4)[[1]]
	#chr,pos1,rsID,refallele,sample,hap,ref,alt
    ii.full = vcf[,c(1:2,i: pmin((i+3), ncol(vcf)))]
    #ii.full<- as.data.frame(vcf[,..cols_keep])
    ii <- ii.full[,-3]
    id=data.frame(ID=unique(ii.full[,3]))
    #nms <- c("CHR", "POS", "HAP", "REF.READS", "ALT.READS")   # Vector of columns you want in this data.frame
    colnames(ii)<-c("CHR", "POS", "HAP", "REF.READS", "ALT.READS")
	#Missing <- setdiff(nms, names(ii))  # Find names of missing columns
	#ii[Missing] <- NA                    # Add them, filled with '0's
	#ii <- ii[nms]                       # Put columns in desired order
	ii.2 = ii[ (ii$HAP == "1|0" | ii$HAP == "0|1") & ii$REF.READS >= opt$min_cov & ii$ALT.READS >= opt$min_cov ,]
    	iii = ii.2[!is.na(ii.2$POS),]# put in phase
al.ref = iii$REF.READS
al.alt = iii$ALT.READS
switch = (iii$HAP == "1|0")
tmp = al.ref[switch]
al.ref[switch] = al.alt[switch]
al.alt[switch] = tmp
if ( !is.na(opt$inp_cnv) ) {
	cnv = read.table( opt$inp_cnv , head=F , as.is=T)
	
	# local parameters to be estimated 
	phi = rep(NA,nrow(cnv))
	mu = rep(NA,nrow(cnv))
	num = rep(NA,nrow(cnv))

	# read through CNVs 
	for ( c in 1:nrow(cnv) ) {
		overlap = iii$CHR == cnv$CHR[c] & iii$POS >= cnv$P0[c] & iii$POS <= cnv$P1[c]
		if ( sum(overlap) > 100 ) {
			fit = vglm(cbind( al.ref[overlap] , al.alt[overlap] ) ~ 1, betabinomialff, trace = FALSE)
			cof = Coef(fit)
			phi[c] = 1/(1+sum(cof))
			mu[c] = cof[1] / sum(cof)
			num[c] = sum(overlap)
		}
	}
	df <- cbind(id,cnv[,c("CHR","P0","P1")],format(cbind(phi,mu),digits=3),num) 
	colnames(df) <- c("ID","CHR","P0","P1","PHI","MU","N")
	write.table(df, quote=F, row.names=F, col.names=!file.exists(paste(opt$out,".local.params",sep='')), append=TRUE, sep='\t' , file=paste(opt$out,".local.params",sep=''))
}

# fit all counts
fit = tryCatch(VGAM::vglm(cbind(al.ref , al.alt) ~ 1, betabinomialff, trace = FALSE), error=function(e) NA)
cof = tryCatch(Coef(fit), error=function(e) data.frame(shape1=NA,shape2=NA))
phi = tryCatch( 1/(1+sum(cof)), error=function(e) NA)
mu = tryCatch(cof[1] / sum(cof), error=function(e) NA)
num = tryCatch(length( al.ref ), error=function(e) NA)
df <- cbind(id,phi , mu, num)
colnames(df) <- c("ID","PHI","MU","N")
write.table(df, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste(opt$out,".global.params",sep='')) , append=TRUE, file=paste(opt$out,".global.params",sep='') )
return(df)
})
} else {
# filter heterozygous sites with minimum reads
vcf = vcf[ (vcf$HAP == "1|0" | vcf$HAP == "0|1") & vcf$REF.READS >= opt$min_cov & vcf$ALT.READS >= opt$min_cov ,]
# put in phase
al.ref = vcf$REF.READS
al.alt = vcf$ALT.READS
switch = (vcf$HAP == "1|0")
tmp = al.ref[switch]
al.ref[switch] = al.alt[switch]
al.alt[switch] = tmp

if ( !is.na(opt$inp_cnv) ) {
	cnv = read.table( opt$inp_cnv , head=F , as.is=T)
	
	# local parameters to be estimated
	phi = rep(NA,nrow(cnv))
	mu = rep(NA,nrow(cnv))
	num = rep(NA,nrow(cnv))

	# read through CNVs 
	for ( c in 1:nrow(cnv) ) {
		overlap = vcf$CHR == cnv$CHR[c] & vcf$POS >= cnv$P0[c] & vcf$POS <= cnv$P1[c]
		if ( sum(overlap) > 100 ) {
			fit = vglm(cbind( al.ref[overlap] , al.alt[overlap] ) ~ 1, betabinomialff, trace = FALSE)
			cof = Coef(fit)
			phi[c] = 1/(1+sum(cof))
			mu[c] = cof[1] / sum(cof)
			num[c] = sum(overlap)
		}
	}
	write.table( cbind(cnv[,c("CHR","P0","P1")],format(cbind(phi,mu),digits=3),num) , quote=F, row.names=F, col.names=c("CHR","P0","P1","PHI","MU","N"), sep='\t' , file=paste(opt$out,".local.params",sep=''))
}

# fit all counts
fit = vglm(cbind( al.ref , al.alt ) ~ 1, betabinomialff, trace = FALSE)              
cof = Coef(fit)
phi = 1/(1+sum(cof))
mu = cof[1] / sum(cof)
num = length( al.ref )
write.table( cbind(phi , mu, num) , quote=F , row.names=F , sep='\t' , col.names=c("PHI","MU","N") , file=paste(opt$out,".global.params",sep='') )
}


