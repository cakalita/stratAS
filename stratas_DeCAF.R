library(VGAM)
library(optparse)
library(plyr)
library(parallel)
library(dplyr)
library(data.table)
setDTthreads(0L)

option_list = list(
	make_option("--input", action="store", default=NA, type='character',
              help="Path to input file [required]"),
	make_option("--samples", action="store", default=NA, type='character',
              help="Path to sample identifier file, must have ID and CONDITION columns [required]"),	
	make_option("--peaks", action="store", default=NA, type='character',
              help="Path to file containing peak/gene boundaries, must contain CHR P0 P1 NAME CENTER columns [required]"),	
	make_option("--global_param", action="store", default=NA, type='character',
              help="Path to global parameter file [required]"),
	make_option("--local_param", action="store", default=NA, type='character',
              help="Path to local parameter file"),
	make_option("--window", action="store", default=100e3 , type='integer',
              help="Window (in bp) for SNPs to test around the peak boundary. [default: %default]"),
	make_option("--perm", action="store", default=0 , type='integer',
              help="# of permutations to shuffle the allele labels (0=off). [default: %default]"),
	make_option("--perm_cond", action="store", default=0 , type='integer',
	            help="# of permutations to shuffle the condition labels (0=off). [default: %default]"),	
	make_option("--min_cov", action="store", default=1 , type='integer',
              help="Individuals must have at least this many reads (for both alleles) to be tested. [default: %default]"),
	make_option("--min_maf", action="store", default=0.01 , type='double',
              help="Minimum minor allele frequency for test SNP. [default: %default]"),
	make_option("--min_het", action="store", default=0.01 , type='double',
              help="Minimum minor heterozygous frequency for test SNP. [default: %default]"),
	make_option("--max_rho", action="store", default=0.10 , type='double',
              help="Maximum local/global over-dispersion parameter for which to include individual in test. [default: %default]"),
	make_option("--cnv_threshold", action="store", default=0.10 , type='double',
              help="Maximum CNV parameter for which to include individual in test. [default: %default]"),
	make_option("--binom", action="store_true", default=FALSE,
              help="Also perform a standard binomial test. [default: %default]"),
	make_option("--bbreg", action="store_true", default=FALSE,
	            help="Also perform a beta binomial regression, requires library(aod). [default: %default]"),	
	make_option("--indiv", action="store_true", default=FALSE,
	            help="Also report the per-individual allele fractions. [default: %default]"),	
	make_option("--exclude", action="store_true", default=75 , type='integer',
              help="The mimium distance between SNPs allowed in the haplotype. [default: %default]"),
	make_option("--gene_express", action="store", default=NA, type='character',
              help="Path to input file containing gene expression results needed for cell-type specific QTL testing [required]"),
	make_option("--sample_name", action="store", default=NA, type='character',
              help="Name of library being analyzed [required]"),
	make_option("--gene_name", action="store", default=NA, type='character',
              help="Name of gene being analyzed [required]"),	
	make_option("--cell_pop", action="store", default=NA, type='character',
              help="Path to file containing cell population fractions, row per individual, columns per cell type [required]"),	
	make_option("--purity", action="store", default=NA, type='character',
              help="Path to file containing tumor purity, row per individual, need column with estimate. ex using Butte et al 2015 CPE column [required]"),	
	make_option("--prs", action="store", default=NA, type='character',
              help="Path to file containing PRS scores, row per individual [required]"),	
	make_option("--out", action="store", default=NA, type='character',
              help="Path to output [default: %default]"),
	make_option("--cell_specific", action="store_true", default=FALSE,
              help="Perform cell-type specific ASE test. [default: %default]"),
	make_option("--data_path", action="store", default=NA, type='character',
              help="Set output data path to store results. [default: %default]"),
	make_option("--eqtl_vanilla", action="store_true", default=FALSE,
	            help="Also perform a standard eqtl test. [default: %default]"),	
	make_option("--bbreg_vanilla", action="store_true", default=FALSE,
	            help="Also perform a beta binomial regression, requires library(aod). [default: %default]"),	
	make_option("--dopurity", action="store_true", default=FALSE,
	            help="Also perform a tumor purity test. [default: %default]"),	
	make_option("--doprs", action="store_true", default=FALSE,
	            help="Also perform a PRS test. [default: %default]"),	
	make_option("--dointeraction", action="store_true", default=FALSE,
	            help="Also perform a eqtl interaction test. [default: %default]"),	
	make_option("--pheno", action="store_true", default=TRUE,
              help="Preforms liklihood ratio test for chisq results between conditions. Curently binary only. [default: %default]"),
	make_option("--ranknorm", action="store_true", default=FALSE,
	            help="rank normalize cell population. [default: %default]"),	
	make_option("--rerun", action="store_true", default=FALSE,
	            help="Restart R script looping from after last snp run if there was a job error. [default: %default]"),		
	make_option("--betabinom", action="store_true", default=FALSE,
              help="Perform a standard beta-binomial test without LRT. [default: %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

PAR.WIN = opt$window
NUM.PERM = opt$perm
NUM.PERM_COND = opt$perm_cond
MIN.MAF = opt$min_maf
MIN.HET = opt$min_het
MIN.COV = opt$min_cov
MAX.RHO = opt$max_rho
DO.BINOM = opt$binom
DO.INDIV = opt$indiv
DO.BBREG = opt$bbreg
DO.CELLSPECIFIC = opt$cell_specific
DO.PHENO = opt$pheno
DO.BETABINOM = opt$betabinom
SAMPLE_NAME = opt$sample_name
GENE_NAME = opt$gene_name
DATA_PATH = opt$data_path
DO.EQTL_VANILLA = opt$eqtl_vanilla
DO.PURITY = opt$dopurity
DO.INTERACTION = opt$dointeraction
DO.PRS = opt$doprs
DO.BBREG_VANILLA = opt$eqtl_vanilla
DO.RERUN = opt$rerun
DO.RANKNORM = opt$ranknorm

message("-----------------------------------")
message("stratAS")
message("https://github.com/gusevlab/stratAS")
message("-----------------------------------")

# --- error checks:
if ( DO.BBREG ) { 
  library("aod")
  if ( is.na(opt$local_param) ) {
    stop("ERROR: --bbreg requires --local_param file\n")
  }
}

if ( NUM.PERM > 0 && NUM.PERM_COND > 0 ) {
  stop("ERROR: --perm and --perm_cond cannot both be set\n")
}

if ( NUM.PERM < 0 || NUM.PERM_COND < 0 ) {
  stop("ERROR: --perm or --perm_cond cannot be negative\n")
}
# ---

message(paste0("Writing to ",DATA_PATH,SAMPLE_NAME))
message("Files being read in")
peaks = read.table( opt$peaks , head=F)
#chr genes peaks file missing header
colnames(peaks) <- c("CHR","P0","P1","NAME","CENTER")
setDF(peaks)
message("peaks:")
head(peaks)
mat = fread( opt$input, fill=T , nThread=1)
setDF(mat)

#the following is just until params
x <- strsplit(SAMPLE_NAME, "[.]")
if (grepl("BRCA",sapply(x, function(y) { y[1] }))) {
samples <- fread(paste0("/agusevlab/ckal/TCGA_vcf/BRCA/samples_BRCA_csv.txt"), header=T, nThread=1)
cnv.all = fread( opt$global_param , head=T, nThread=1)
setDF(cnv.all)
cnv.all <- transform(cnv.all, barcode=paste0(ID,"-",sample_type))
} else {
samples <- fread(paste0("/agusevlab/ckal/TCGA_vcf/KIRC/samples_KIRC_new.txt"), header=T, nThread=1)
samples = unique(subset(samples, sample_type=="01A"))
cnv.all = fread( opt$global_param , head=T, nThread=1)
setDF(cnv.all)
cnv.all = unique(subset(cnv.all, sample_type=="01A"))
cnv.all <- transform(cnv.all, barcode=paste0(ID,"-",sample_type))
}
message("samples:")
head(samples)
#adds individual ids to genotype/counts columns in the ASE matrix
samples_rep <- lapply(split(samples,rownames(samples)),function(i){
	#i <- split(samples,rownames(samples))[[1]]
	ii<-rep(paste0(i$ID,"-",i$sample_type),4)
	#iii <- as.vector(ldply(ii,data.frame)[,2])
	return(ii)})
samples_repv <- unlist(samples_rep,use.names = FALSE)

samples_repv_mat <- c("CHR","POS","RSID","REF","ALT",samples_repv)
colnames(mat) <- samples_repv_mat

if (grepl("BRCA",sapply(x, function(y) { y[1] }))) {
mat <- mat[!grepl("-11A|-01B|-06A|-11B",colnames(mat))]
samples_repv_s <- samples_repv_mat[!grepl("-11A|-01B|-06A|-11B",samples_repv_mat)]
samples_repv_mat <- samples_repv_s
colnames(mat) <- samples_repv_mat
}

message("matrix:")
head(mat[,c(1:10)])

#########
cnv.all <- cnv.all[cnv.all$barcode %in% samples_repv_mat,]
message("global cnv:")
head(cnv.all)

if (DO.PHENO) {
	phe = read.table( opt$samples , head=T , as.is=T)
	message("Doing condition test")
}
if (DO.CELLSPECIFIC | DO.PURITY | DO.PRS) {
	gene_express = fread( opt$gene_express , head=T , sep='\t', check.names=FALSE, nThread=1) #tcga samples in header were getting converted from dast to dot
	setDF(gene_express)
	gene_express <- gene_express[colnames(gene_express) %in% c(unique(samples_repv_mat),"NAME")]
	message("gene express:")
	head(gene_express)
	cell_pop.both = fread( opt$cell_pop , head=T , nThread=1)
	setDF(cell_pop.both)
	if (grepl("KIRC|BRCA",sapply(x, function(y) { y[1] }))) {
	cell_pop = unique(subset(cell_pop.both, sample_type=="01A"))
	} else {
	cell_pop = cell_pop.both
	}
	cell_pop <- cell_pop[paste0(cell_pop$sample_ID,"-",cell_pop$sample_type) %in% samples_repv_mat,]
	rm(cell_pop.both)
	message("cell pop:")
	head(cell_pop[,c(1:5)])
	x <- strsplit(SAMPLE_NAME, "[.]")
	if(grepl("xcell",opt$cell_pop)){
		message("running xcell")
		xcell_midcommon_file <- fread(file="/agusevlab/ckal/TCGA_vcf/BRCA/xcell_midcommon_BRCA.txt",header=F, nThread=1)
	}
	message("Doing condition test")
}
if(DO.PURITY | DO.CELLSPECIFIC ) {
	purity_in = fread( opt$purity , head=T , sep='\t', nThread=1) 
	purity_df1 = purity_in[,c("sample_ID","CPE")]
	setDF(purity_df1)
	purity_df <- purity_df1[purity_df1$sample_ID %in% samples_repv_mat,]
	rm(purity_in)
	message("purity:")
	head(purity_df)
	message("Doing tumor purity test")
}
if (DO.PRS) {
#/agusevlab/ckal/PRS/PRC_KIRC.txt
	PRC <- fread(opt$prs, nThread=1)
	message("Doing PRS test")
}

message("Files have been read in")

PERM.PVTHRESH = 0.05 / nrow(mat)

bb.loglike = function( mu , ref , alt, rho ) {
        keep = !is.na(ref + alt) & ref + alt > 0
        -1 * sum( dbetabinom(alt[keep], (ref+alt)[keep], mu, rho = rho[keep], log = T) )
}

bbinom.test = function( ref , alt , rho ) {
	if ( length(ref) > 0 && length(alt) > 0 && length(rho) > 0 ) {
		opt = optimize( bb.loglike , interval=c(0,1) , ref , alt, rho )
		opt$lrt = 2 * (opt$objective - bb.loglike( 0.5 , ref , alt , rho) )
		opt$pv = pchisq( abs(opt$lrt) , df=1 , lower.tail=F )
	} else {
		opt = list( "lrt" = NA , "pv" = NA , "min" = NA , "objective" = NA )
	}
	return( opt )
}

bbreg.cf0.test = function( ref , alt , rho , purity=NULL, cond=NULL ) {
    if ( !is.null(cond) & !is.null(purity)) {
      df = data.frame( y=alt , n=ref+alt , cond = cond , purity=purity, rho=rho )
    } else if ( !is.null(cond) & is.null(purity)) {
      df = data.frame( y=alt , n=ref+alt , cond = cond , rho=rho )  
    } else {
      df = data.frame( y=alt , n=ref+alt , rho=rho )
    }
    # test with a fixed overd parameter for each individual (for some reason this is very SLOW!)
    # reg = betabin( cbind( y , n - y  ) ~ 1 + cond , ~ phi.group , df , fixpar = list( 3:(n+2) , rho ) )
    df = df[ df$n > 0, ]
    
    # for efficiency, discretize the overd parameters into five groups
    nq = min( 5 , length(unique(rho)) )
    phi.q = rep(NA,nrow(df))
    rho.q = rep(NA,nq)
    qq = quantile(unique(df$rho), probs = seq(0, 1, 1/nq))
    for ( i in 1:nq ) {
      keep = df$rho >= qq[i] & df$rho <= qq[i+1]
      phi.q[ keep ] = i
      rho.q[ i ] = mean( df$rho[keep] )
    }
    df$phi.q = as.factor(phi.q)

    statement= nrow(df) < 2 || sd(df$cond,na.rm=T) == 0
    # need at least two samples to do regression with covariates
    if ( !is.null(cond) & !is.null(purity)) {
      # need at least two samples in both conditions
      if ( statement) {
        return( list( "pv" = c(NA,NA,NA) ) ) 
        break
    } else reg = betabin( cbind( y , n-y  ) ~ 1 + cond + purity , ~ phi.q , df , fixpar = list( 4:(nq+3) , rho.q ) )
    } else if (!is.null(cond) & is.null(purity)) {
      if ( statement) {
        return( list( "pv" = c(NA,NA) ) ) 
        break
    } else reg = betabin( cbind( y , n - y  ) ~ 1 + cond , ~ phi.q , df , fixpar = list( 3:(nq+2) , rho.q ) )
    } else {
        reg = betabin( cbind( y, n-y ) ~ 1 , ~ phi.q , df , fixpar = list( 2:(nq+1) , rho.q ) )
    }
    
    zscores = coef(reg) / sqrt(diag(vcov( reg )))
    pvals = 2*(pnorm( abs(zscores) , lower.tail=F))
    opt = list( "pv" = pvals, "z" = zscores )
  return( opt )
}

rank.normalize <- function(x, FUN=qnorm, ties.method = "average", na.action) {
	if (missing(na.action)) {
		na.action <- get(getOption("na.action"))
	}
	if(! is.function(na.action)) {
		stop("'na.action' must be a function")
	}
	x <- na.action(x)
	FUN(rank(x, ties.method = ties.method)/(length(x)+1))
}

cur.chr = unique(mat[,1])
m = match(peaks$CHR , cur.chr )
peaks = peaks[!is.na(m),]

message(paste0(opt$gene_name))
if ( !is.null(opt$gene_name) &!is.na(opt$gene_name)) {
	message(paste0(opt$gene_name, " is not null"))
	peaks = na.omit(peaks[peaks$NAME==GENE_NAME,])
	SAMPLE_NAME.x <- strsplit(SAMPLE_NAME, "[.]")
	SAMPLE_NAME.old <- SAMPLE_NAME
	SAMPLE_NAME <- paste(paste0(sapply(SAMPLE_NAME.x, function(y){y[1:4]}),collapse="."),GENE_NAME,paste0(sapply(SAMPLE_NAME.x, function(y){y[5:6]}),collapse="."),sep=".")
	message(paste0("peaks  = ",peaks$NAME))
} else {
	peaks = na.omit(peaks)
head(peaks)
}

if ( !is.na(opt$local_param) )  {
	cnv.local = read.table( opt$local_param , head=T ,as.is=T,fill=T)
	cnv.local <- cnv.local[cnv.local$ID %in% c(cnv.all$ID),]
	m_cnv = match(cnv.local$CHR , cur.chr )
	if (grepl("KIRC|BRCA",sapply(x, function(y) { y[1] }))) {
	cnv.local = cnv.local[!is.na(m_cnv) & cnv.local$sample_type=="01A",]
	} else {
	cnv.local = cnv.local[!is.na(m_cnv),]
	}	
	cnv.local <- cnv.local[!duplicated(cnv.local[,c("ID","CHR","P0","P1")]),]
	cnv.local <- transform(cnv.local, barcode=paste0(ID,"-",sample_type))
	rm(m_cnv)
	message("cnvlocal:")
	head(cnv.local)
}

message("starting matrix")

N = (ncol(mat) - 5)/4
M = nrow(mat)

HAPS = list()
GENO.H1 = matrix( 0 , nrow=M , ncol=N )
GENO.H2 = matrix( 0 , nrow=M , ncol=N )

if (DO.PHENO) {
	PHENO = phe$CONDITION
	m_phe = match( phe$ID , cnv.all$ID )
	cnv.all = cnv.all[m_phe,]
	RHO.ALL = cnv.all$PHI
	rm(m_phe)
} else {
	RHO.ALL = cnv.all$PHI
}

message("starting haps")
for ( h in 1:2 ) {
	HAPS[[h]] = matrix( 0 , nrow=M , ncol=N )
}

#head(mat[,c(1:10)])
message("starting genos")
# standardize matrix to the same haplotypes
for ( i in 1:N ) {
	GENO.H1[,i] = mat[ , 6 + 4*(i-1) ]
	GENO.H2[,i] = mat[ , 6 + 4*(i-1) + 1 ]
	HET = GENO.H1[,i] != GENO.H2[,i]
	cur.ALT = mat[ , 6 + 4*(i-1) ] == 0
	HAPS[[1]][HET & cur.ALT,i] = mat[ HET & cur.ALT , 6 + 4*(i-1) + 2 ]
	HAPS[[1]][HET & !cur.ALT,i] = mat[ HET & !cur.ALT , 6 + 4*(i-1) + 3 ]
	HAPS[[2]][HET & cur.ALT,i] = mat[ HET & cur.ALT , 6 + 4*(i-1) + 3 ]
	HAPS[[2]][HET & !cur.ALT,i] = mat[ HET & !cur.ALT , 6 + 4*(i-1) + 2 ]	
}

message("matrix made")

#save(GENO.H1,GENO.H2, file=paste0(DATA_PATH,SAMPLE_NAME,"GENOobjects.RData"))
#saveRDS(HAPS, file=paste0(DATA_PATH,SAMPLE_NAME,".HAPS.rds") , compress=FALSE)

#rm(HAPS,GENO.H1,GENO.H2)

options( digits = 4 )

RHO = RHO.ALL
RHO[ RHO > MAX.RHO ] = NA

rm(RHO.ALL,M,cur.chr)
gc()

COL.HEADER = c("CHR","POS","RSID","P0","P1","NAME","CENTER","N.HET","N.READS")
COL.HEADER.PHENO = c("ALL.AF","ALL.BBINOM.P","C0.AF","C0.BBINOM.P","C1.AF","C1.BBINOM.P","DIFF.BBINOM.P")
COL.HEADER.BINOM = c("ALL.BINOM.P","ALL.C0.BINOM.P","ALL.C1.BINOM.P","FISHER.OR","FISHER.DIFF.P")
COL.HEADER.PHENO.INDIV = c("IND.C0","IND.C0.COUNT.REF","IND.C0.COUNT.ALT","IND.C1","IND.C1.COUNT.REF","IND.C1.COUNT.ALT")
COL.HEADER.INDIV = c("IND","IND.COUNT.REF","IND.COUNT.ALT")
COL.HEADER.BBREG = c("C0.BBREG.P","C0.CNV.BBREG.P","C1.BBREG.P","C1.CNV.BBREG.P","ALL.BBREG.P","DIFF.BBREG.P","DIFF.CNV.BBREG.P")
COL.HEADER.BETABINOM = c("ALL.AF","ALL.BBINOM.P")

HEAD = COL.HEADER

if (DO.RERUN & !DO.CELLSPECIFIC & !DO.PRS) {
	cfqtl_file <- read.table(file=paste0(DATA_PATH,SAMPLE_NAME,".eqtl_vanilla_stratas_results.txt"), header=T, fill=T)
	cfqtl <- subset(cfqtl_file, NAME %in% peaks$NAME)
	last_snp <- tail(cfqtl,2)
	last_snp <- subset(last_snp, !RSID=="RSID")
if (length(last_snp$RSID)==2) {
	last_snp <- last_snp[2,]
}
 	peaks= peaks[which(peaks$CHR==last_snp$CHR & peaks$P0==last_snp$P0 & peaks$P1==last_snp$P1)+1:nrow(peaks), ]
}

if (DO.RERUN & DO.CELLSPECIFIC) {
	cfqtl_file <- read.table(file=paste0(DATA_PATH,SAMPLE_NAME,".cfQTL.bbreg_stratas_results.txt"), header=T, fill=T)
	cfqtl <- subset(cfqtl_file, NAME %in% peaks$NAME)
	last_snp <- tail(cfqtl,2)
	last_snp <- subset(last_snp, !RSID=="RSID")
if (length(last_snp$RSID)==2) {
	last_snp <- last_snp[2,]
}
 	peaks= peaks[which(peaks$CHR==last_snp$CHR & peaks$P0==last_snp$P0 & peaks$P1==last_snp$P1)+1:nrow(peaks), ]
}
if (DO.RERUN & DO.PRS) {
	prsfile1 <- fread(file=paste0(DATA_PATH,SAMPLE_NAME,".PRSQTL_stratas_results.txt"), nThread=1)
	prsfile <- subset(prsfile1, NAME %in% peaks$NAME)
	last_snp <- tail(prsfile,2)
	last_snp <- subset(last_snp, !RSID=="RSID")
if (length(last_snp$RSID)==2) {
	last_snp <- last_snp[2,]
}
 	peaks= peaks[which(peaks$CHR==last_snp$CHR & peaks$P0==last_snp$P0 & peaks$P1==last_snp$P1)+1:nrow(peaks), ]
}

for ( p in 1:nrow(peaks) ) {
	#message(paste0("Running ", p))
	#for testing uncomment out next line
	#p=2
	#browser()
	if (DO.CELLSPECIFIC | DO.PURITY | DO.PRS | DO.EQTL_VANILLA) {
		gene = na.omit(merge(peaks[p,], gene_express, by="NAME"))
		message(paste0("gene  = ",gene$NAME))
		Y.total = gene[,-c(1:5)]
		message(paste0("Running ", p, "=", gene[,"NAME"]))
	}

	if ( DO.RERUN ) {
		mat = mat[which(mat$V3==last_snp$RSID)+1:nrow(mat),]
	}

	cur = mat[,1] == peaks$CHR[p] & mat[,2] >= peaks$P0[p] & mat[,2] <= peaks$P1[p]
	cur[is.na(cur)] <- FALSE

	if(sum(cur) == 0 ) next()

	if ( !is.na(opt$local_param) )  {
		cnv.a <- merge(cnv.all,cnv.local, by="barcode", all=T)
		cnv.a <- transform(cnv.a, PHI=ifelse(is.na(PHI.y), PHI.x, PHI.y))
		cur.cnv_p = cnv.a[ cnv.a$CHR == peaks$CHR[p] & cnv.a$P0 < peaks$P0[p] & cnv.a$P1 > peaks$P1[p] , ]
		cur.cnv <- merge(cur.cnv_p, data.frame(ID=cnv.all$ID, sample_type=cnv.all$sample_type, barcode=cnv.all$barcode),all.y=T)
	
		if (DO.PHENO) {
			m = match( phe$ID , cur.cnv$ID )
			cur.cnv = cur.cnv[m,]
		}
		cur.cnv$PHI <- ifelse(cur.cnv$PHI>MAX.RHO, NA, cur.cnv$PHI)
		COVAR = cur.cnv$CNV
		COVAR[is.na(COVAR)] <- mean(COVAR,na.rm=T)
		COVAR[ COVAR > as.numeric(opt$cnv_threshold) ] = NA
		cur.cnv$CNV=COVAR
		cur.cnv.l <- ldply(lapply(split(cur.cnv,cur.cnv$barcode), function(i){
			#i=split(cur.cnv,cur.cnv$ID)[[99]] i=subset(cur.cnv,ID=="TCGA-BH-A0HO")
			if(all(is.na(i$CNV))){
				bestCNV="NA"
				ii <- i[1,]
				#message(paste0(i$barcode," all cnv NAs"))
			} else{
			bestCNV=min(abs(i$CNV),na.rm=T)
			ii <- i[abs(i$CNV)==bestCNV,]
			ii<-ii[!duplicated(ii$barcode),]
			#message(paste0(i$barcode," taking min CNV"))
			}
			return(ii[!is.na(ii$barcode),])
			}),data.frame)
		RHO = cur.cnv.l$PHI
		COVAR = cur.cnv.l$CNV
		IDS = subset(cur.cnv.l, !is.na(CNV))$barcode
	} else {
		IDS = cnv.all[,"barcode"]
	}

	if (DO.CELLSPECIFIC | DO.PURITY | DO.PRS | DO.EQTL_VANILLA) {
		cnv.ids <- cnv.all[,c("ID","sample_type","barcode")]
		cell_pop <- transform(cell_pop, barcode=paste0(sample_ID,"-",sample_type))
		cell_pop_sub1 <- merge(unique(cnv.ids[,c("sample_type","barcode")]), cell_pop, by=c("sample_type","barcode"),all.x=T)
		eqtl_IDset <- ifelse(unique(cell_pop_sub1$barcode) %in% colnames(Y.total) & unique(cell_pop_sub1$barcode) %in% IDS, TRUE, FALSE)
		cell_pop_sub <- subset(cell_pop_sub1, barcode %in% colnames(Y.total) & barcode %in% IDS)
		cell_pop_sub <- cell_pop_sub[!duplicated(cell_pop_sub$barcode),]
		cur.y.total = colnames(Y.total) %in% IDS
		Y.total.ind <- unname(unlist(Y.total[cur.y.total]))
		Y.total.ind.norm = rank.normalize(Y.total.ind)
		cell_types <- colnames(cell_pop_sub[,-c(1:3)])
		if(grepl("xcell",opt$cell_pop)){
			cell_types <- cell_types[cell_types %in% xcell_midcommon_file$V1]
		}
	}
	if (DO.PURITY | DO.CELLSPECIFIC ) {
		purity_df = transform(purity_df, NF=1-CPE)
		purity_sub2 <- merge(purity_df, cnv.ids, by.x="sample_ID",by.y="barcode",all.y=T)
		purity_sub1 <- purity_sub2[!duplicated(purity_sub2$sample_ID),]
		purity_sub <- subset(purity_sub1, sample_ID %in% names(unlist(Y.total[cur.y.total])) & sample_ID %in% IDS)
		purityIDS <- IDS[IDS %in% purity_sub$sample_ID]
		cur.y.total = colnames(Y.total) %in% purityIDS
		Y.total.ind <- unname(unlist(Y.total[cur.y.total]))
		Y.total.pur.norm = rank.normalize(Y.total.ind)
		eqtlpur_IDset <- ifelse(unique(cell_pop_sub1$barcode) %in% purityIDS, TRUE, FALSE)
		NF_qtl = purity_sub$NF
		NF_qtl[is.na(NF_qtl)] <-mean(NF_qtl,na.rm=T)
	}
	if (DO.PRS) {
		#this needs updating
		PRS_sub1 <- merge(cnv.ids, PRC, by.x="ID", by.y="sample_ID",all.x=T)
		PRS_sub1 <- subset(PRS_sub1, ID %in% colnames(Y.total))
		PRS_sub <- PRS_sub1[!duplicated(PRS_sub1$ID),]
		PRS_score = PRS_sub$SCORE1_AVG
		PRS_score[is.na(PRS_score)] <-mean(PRS_score,na.rm=T)
	}

	# collapse reads at this peak
	cur.h1 = vector()
	cur.h2 = vector()
	cur.i = vector()

	#load(file=paste0(DATA_PATH,SAMPLE_NAME,"GENOobjects.RData"))
	#HAPS <- readRDS(file=paste0(DATA_PATH,SAMPLE_NAME,".HAPS.rds"))

	for ( i in 1:N ) {
		reads.keep = GENO.H1[cur,i] != GENO.H2[cur,i] & HAPS[[1]][cur,i] >= MIN.COV & HAPS[[2]][cur,i] >= MIN.COV
		reads.keep[reads.keep] <- !c(FALSE, diff(mat[cur, 2][reads.keep]) < opt$exclude)
		cur.h1 = c( cur.h1 , (HAPS[[1]][cur,i])[reads.keep] )
		cur.h2 = c( cur.h2 , (HAPS[[2]][cur,i])[reads.keep] )
		cur.i = c( cur.i , rep( i , sum(reads.keep)) )
	}

	#rm(HAPS,GENO.H1,GENO.H2)

	if ( length(unique(cur.i)) > MIN.MAF*N && sum(cur.h1) + sum(cur.h2) < 0 ) {
		message("length(unique(cur.i)) > MIN.MAF*N && sum(cur.h1) + sum(cur.h2) < 0")
	} else {
		# test all nearby SNPs

		if( PAR.WIN == -1 ) {
			cur.snp = mat[,1] == peaks$CHR[p] & mat[,2] >= peaks$P0[p] & mat[,2] <= peaks$P1[p]
		} else {
			cur.snp = mat[,1] == peaks$CHR[p] & mat[,2] >= peaks$CENTER[p] - PAR.WIN & mat[,2] <= peaks$CENTER[p] + PAR.WIN
		}

		snp_list <- lapply(which(cur.snp), function(s){
	#for testing uncomment out next line
		#s = 1

			#load(file=paste0(DATA_PATH,SAMPLE_NAME,"GENOobjects.RData"))
			
			mat_s = mat[s,]
			GENO.H1_s = GENO.H1[s,]
			GENO.H2_s = GENO.H2[s,]
			#rm(GENO.H1,GENO.H2)

			# restrict to hets
			HET = GENO.H1_s != GENO.H2_s & !is.na(RHO) 
			DO.HET = sum(HET) > MIN.HET*N
			message(paste0("DO.HET = ",DO.HET))
	
			# collect REF/ALT heterozygous haplotypes
				m1 = match( cur.i , which(HET & GENO.H1_s == 0 & !is.na(COVAR)) )
				m2 = match( cur.i , which(HET & GENO.H2_s == 0 & !is.na(COVAR)) )
				CUR.REF = c( cur.h1[ !is.na(m1) ] , cur.h2[ !is.na(m2) ])
				CUR.ALT = c( cur.h2[ !is.na(m1) ] , cur.h1[ !is.na(m2) ])
				DO.REFALT = sum( (CUR.REF+CUR.ALT) > 0 , na.rm=T) > 1
				message(paste0("DO.REFALT = ",DO.REFALT))

				CUR.IND = c( cur.i[ !is.na(m1) ] , cur.i[ !is.na(m2)] )
				CUR.ID = unique(IDS)[CUR.IND] #this is on the theory that CUR.IND is vector of Nth individuals representing that snp
				df_ind <- cbind(paste(CUR.IND,collapse=',') , paste(CUR.REF,collapse=',') , paste(CUR.ALT,collapse=',') , paste(CUR.ID,collapse=','))
				df_info <- cbind(mat_s[,1:3] , peaks[p,c("P0","P1","NAME","CENTER")] , sum(HET) , sum(CUR.REF) + sum(CUR.ALT))
		          	colnames(df_info) <- COL.HEADER	
				all.RHO = RHO[CUR.IND]
				if (DO.PURITY | DO.CELLSPECIFIC ) {
				NF_ase = purity_sub1$NF
				NF_ase[is.na(NF_ase)] <-mean(NF_ase,na.rm=T)
  				all.NF_ase = NF_ase[CUR.IND]
				all.NF_ase[is.na(all.NF_ase)] <-mean(all.NF_ase,na.rm=T)
				}
				if (DO.PRS) {
  				all.PRS_score = PRS_score[CUR.IND]
				all.PRS_score[is.na(all.PRS_score)] <-mean(all.PRS_score,na.rm=T)
				}

				rm(HET)
				gc()

				if (DO.EQTL_VANILLA | DO.PURITY | DO.INTERACTION | DO.PRS) {
				SNP_geno=GENO.H1_s + GENO.H2_s
				SNP <- SNP_geno[eqtl_IDset]
				SNPpur <- SNP_geno[eqtlpur_IDset]
				}
				#purity
				if (DO.PURITY & DO.INTERACTION) {
				message("Doing Purity QTL test")
				lm_eqtl_purity = tryCatch(lm(Y.total.pur.norm ~ NF_qtl*SNPpur + SNPpur + NF_qtl ), error=function(e) NA)
       				lm_eqtl_purity_df = tryCatch(as.data.frame(summary(lm_eqtl_purity)$coefficients), error=function(e) NA)
   				lm_eqtl_purity_interaction = tryCatch(lm_eqtl_purity_df["SNPpur",c(3:4)], error=function(e) data.frame(z_eqtl=NA,z_eqtl_pval=NA))
    				colnames(lm_eqtl_purity_interaction) <- c("z_eqtl","z_eqtl_pval")
    				z_eqtl_purity = tryCatch(data.frame(z_eqtl=lm_eqtl_purity_interaction$z_eqtl, z_eqtl_pval=2*pnorm(-abs(lm_eqtl_purity_interaction$z_eqtl))), error=function(e) data.frame(z_eqtl=NA,z_eqtl_pval=NA))
    				purity_qtl <- tryCatch(cbind(df_info, z_eqtl_purity , SAMPLE_NAME), error=function(e) data.frame(df_info,beta=NA,se=NA,z_eqtl=NA,z_eqtl_pval=NA,SAMPLE_NAME))
				purity_qtl <- purity_qtl %>% mutate_all(na_if,"")
				rm(lm_eqtl_purity,lm_eqtl_purity_df,lm_eqtl_purity_interaction,z_eqtl_purity)
				}

				#PRS
				if (DO.PRS) {
				message("Doing PRS QTL test")
				lm_eqtl_PRS = tryCatch(lm(Y.total.ind.norm ~ PRS_score*SNP + SNP + PRS_score + NF_qtl), error=function(e) NA)
       				lm_eqtl_PRS_df = tryCatch(as.data.frame(summary(lm_eqtl_PRS)$coefficients), error=function(e) NA)
   				lm_eqtl_PRS_interaction = tryCatch(lm_eqtl_PRS_df["PRS_score:SNP",], error=function(e) data.frame(beta=NA,se=NA,z_eqtl=NA,z_eqtl_pval=NA))
    				colnames(lm_eqtl_PRS_interaction) <- c("beta","se","z_eqtl","z_eqtl_pval")
    				z_eqtl_PRS = tryCatch(data.frame(beta=lm_eqtl_PRS_interaction$beta,se=lm_eqtl_PRS_interaction$se, z_eqtl=lm_eqtl_PRS_interaction$z_eqtl, z_eqtl_pval=2*pnorm(-abs(lm_eqtl_PRS_interaction$z_eqtl))), error=function(e) data.frame(beta=NA,se=NA,z_eqtl=NA,z_eqtl_pval=NA))
    				PRS_qtl <- tryCatch(cbind(df_info, z_eqtl_PRS , SAMPLE_NAME), error=function(e) data.frame(df_info,beta=NA,se=NA,z_eqtl=NA,z_eqtl_pval=NA,SAMPLE_NAME))
				PRS_qtl <- PRS_qtl %>% mutate_all(na_if,"")
				rm(lm_eqtl_PRS,lm_eqtl_PRS_df,lm_eqtl_PRS_interaction,z_eqtl_PRS )
				}

				if (DO.EQTL_VANILLA) {
				if (DO.PURITY) {
				lm_eqtl_vanilla = tryCatch(lm(Y.total.pur.norm ~ SNPpur + NF_qtl ), error=function(e) NA)
				lm_eqtl_vanilla_df = tryCatch(as.data.frame(summary(lm_eqtl_vanilla)$coefficients), error=function(e) NA)
   				lm_eqtl_vanilla_interaction = tryCatch(lm_eqtl_vanilla_df["SNPpur",c(3,4)], error=function(e) data.frame(z_eqtl_vanilla=NA,eqtl_pval_vanilla=NA))
       				} else {
				lm_eqtl_vanilla = tryCatch(lm(Y.total.ind.norm ~ SNP ), error=function(e) NA)
				lm_eqtl_vanilla_df = tryCatch(as.data.frame(summary(lm_eqtl_vanilla)$coefficients), error=function(e) NA)
   				lm_eqtl_vanilla_interaction = tryCatch(lm_eqtl_vanilla_df["SNP",c(3,4)], error=function(e) data.frame(z_eqtl_vanilla=NA,eqtl_pval_vanilla=NA))
				}
    				colnames(lm_eqtl_vanilla_interaction) <- c("z_eqtl_vanilla","eqtl_pval_vanilla")
    				z_eqtl_vanilla = tryCatch(data.frame(z_eqtl_vanilla=lm_eqtl_vanilla_interaction$z_eqtl_vanilla,eqtl_pval_vanilla=2*pnorm(-abs(lm_eqtl_vanilla_interaction$z_eqtl_vanilla))), error=function(e) data.frame(z_eqtl_vanilla=NA,eqtl_pval_vanilla=NA))
				fwrite(cbind(df_info,z_eqtl_vanilla, SAMPLE_NAME), quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".eqtl_vanilla_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".eqtl_vanilla_stratas_results.txt") )
				#message(paste0(lm_eqtl_vanilla_interaction$z_eqtl_vanilla,": test finished"))
				rm(lm_eqtl_vanilla,lm_eqtl_vanilla_df,lm_eqtl_vanilla_interaction,z_eqtl_vanilla)
				}				

				if (DO.CELLSPECIFIC) {
				message("Doing Interaction QTL test")
				eqtl_list <- mclapply(cell_types, function(c){
					#c=cell_types[1]
       				cf0_cell <- cell_pop_sub[, c]
       				cf0 <- cf0_cell[eqtlpur_IDset]
					cf0[is.na(cf0)] <-mean(cf0,na.rm=T)
					if (DO.RANKNORM) {
						cf0=rank.normalize(cf0)
					}
 					lm_eqtl = tryCatch(lm(Y.total.pur.norm ~ cf0*SNPpur + SNPpur + cf0 + NF_qtl + NF_qtl*SNPpur), error=function(e) NA)
       					lm_eqtl_df = tryCatch(as.data.frame(summary(lm_eqtl)$coefficients), error=function(e) NA)
   					lm_interaction = tryCatch(lm_eqtl_df["cf0:SNP",c(3,4)], error=function(e) data.frame(z_eqtl=NA,pval=NA))
    					colnames(lm_interaction) <- c("z_eqtl","pval")
    					z_eqtl = tryCatch(data.frame(cell=c, z_eqtl=lm_interaction$z_eqtl,z_eqtl_pval=2*pnorm(-abs(lm_interaction$z_eqtl))), error=function(e) data.frame(cell=c,z_eqtl=NA,z_eqtl_pval=NA))
					df <- cbind(df_info, z_eqtl)
					return(df)})
				z_eqtl_df <- ldply(eqtl_list, data.frame)
				rm(SNP,SNP_geno,eqtl_list)
				#message(paste0(s,": eqtl test finished"))
 				}	

			if ( DO.HET & DO.REFALT ) {
				if ( !sum(CUR.REF) + sum(CUR.ALT) > 1 ) {
					message("!sum(CUR.REF) + sum(CUR.ALT) > 1")
				} else {

					for ( perm in 0:max(NUM.PERM,NUM.PERM_COND) ) {
					  if (DO.PHENO) {
					  	CUR.PHENO = PHENO
					  	}
					  
						if ( perm > 0 ) {
							# randomly swap REF/ALT alleles
						  if ( NUM.PERM > 0 ) {
  							cur.swap = unique(CUR.IND)[ as.logical(rbinom( length(unique(CUR.IND)) , 1 , 0.5 )) ]
  							cur.swap = !is.na(match( CUR.IND , cur.swap ))
  							tmp = CUR.REF[ cur.swap ]
  							CUR.REF[ cur.swap ] = CUR.ALT[ cur.swap ]
  							CUR.ALT[ cur.swap ] = tmp
  						# shuffle the phenotype label
  						} else if (NUM.PERM_COND > 0 ) {
  							if (DO.PHENO) {
  								CUR.PHENO = sample(PHENO)
  							}
						  }
						}
						if ( DO.BETABINOM ) {
						# --- preform beta-binomial test without LRT to compare conditions
							HEAD = c(COL.HEADER,COL.HEADER.BETABINOM,"SAMPLE_NAME")
							tst.bbinom.ALL = tryCatch(bbinom.test( CUR.REF , CUR.ALT , RHO[CUR.IND] ), error=function(e) NA)
							df <- tryCatch(cbind(df_info , tst.bbinom.ALL$min , tst.bbinom.ALL$pv , SAMPLE_NAME), error=function(e) data.frame(df_info,AF=tst.bbinom.ALL$min,pval=tst.bbinom.ALL$pv,SAMPLE_NAME))
							colnames(df) <- HEAD
							fwrite(df, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".bbinom_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".bbinom_stratas_results.txt") )
							#message(paste0(s,": beta.binom test finished"))
							rm(tst.bbinom.ALL,df)
						# --- print individual counts ##the testing loop needs to be fixed
							if ( DO.INDIV ) {
								test_df <- data.frame(ID=CUR.ID,ref=CUR.REF,alt=CUR.ALT,rho=all.RHO)
								test_list <- ldply(lapply(split(test_df,rownames(test_df)), function(i){
  									#i<-split(test_df,rownames(test_df))[[1]]
  									tst.bbinom.ALL = tryCatch(bbinom.test( i$ref , i$alt , i$rho ), error=function(e) NA)
									df <- tryCatch(cbind(df_info , ID=i$ID,ref=i$ref,alt=i$alt, AF=tst.bbinom.ALL$min , pval=tst.bbinom.ALL$pv , SAMPLE_NAME), error=function(e) data.frame(df_info,ID=i$ID,ref=i$ref,alt=i$alt,AF=NA,pval=NA,SAMPLE_NAME))
									df <- transform(df, beta=AF-0.5, z=qnorm(pval))
									df <- transform(df, z=ifelse(beta<0, abs(z)*-1, abs(z)))
									df <- transform(df, se=beta/z)
									return(df)
  								}),data.frame)
								fwrite(test_list, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".bbinom.ind_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".bbinom.ind_stratas_results.txt") )
							}
						}
						if ( DO.BBREG_VANILLA & !DO.PHENO ) {
							tst.bbreg.vanilla = tryCatch(bbreg.cf0.test( CUR.REF , CUR.ALT , all.RHO ), error=function(e) NA)
							tst.bbreg.vanilla_res = tryCatch(as.data.frame(tst.bbreg.vanilla)["(Intercept)",c(2,1)], error=function(e) data.frame(z_AI=NA,z_AI_pval=NA))
							colnames(tst.bbreg.vanilla_res) <- c("z_AI","z_AI_pval")		
							df <- tryCatch(cbind(df_info, tst.bbreg.vanilla_res , SAMPLE_NAME), error=function(e) data.frame(df_info,z_AI=NA,z_AI_pval=NA,SAMPLE_NAME))
							fwrite(df, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".bbreg_vanilla_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".bbreg_vanilla_stratas_results.txt") )
							#message(paste0(s,": bbreg vanilla test finished"))
							rm(tst.bbreg.vanilla,tst.bbreg.vanilla_res,df)
						}
						if (DO.BBREG & DO.PRS ) {
							tst.bbreg.PRS = tryCatch(bbreg.cf0.test( CUR.REF , CUR.ALT , rho=all.RHO , purity=all.NF_ase, cond=all.PRS_score), error=function(e) NA)
							tst.bbreg.PRS_res = tryCatch(as.data.frame(tst.bbreg.PRS)["cond",c(2,1)], error=function(e) data.frame(z_AI=NA,z_AI_pval=NA))
							colnames(tst.bbreg.PRS_res) <- c("z_AI","z_AI_pval")
							PRS_ase <- tryCatch(cbind(df_info[c("RSID")], tst.bbreg.PRS_res), error=function(e) data.frame(df_info[c("RSID")],z_AI=NA,z_AI_pval=NA))
							#setDF(PRS_ase)
							PRS_ase <- PRS_ase %>% mutate_all(na_if,"")
							rm(tst.bbreg.PRS,tst.bbreg.PRS_res)
						}
						if (DO.PURITY & DO.INTERACTION) {
							tst.bbreg.purity = tryCatch(bbreg.cf0.test( CUR.REF , CUR.ALT , rho=all.RHO , cond=all.NF_ase), error=function(e) NA)
							tst.bbreg.purity_res = tryCatch(as.data.frame(tst.bbreg.purity)["(Intercept)",c(2,1)], error=function(e) data.frame(z_AI=NA,z_AI_pval=NA))
							colnames(tst.bbreg.purity_res) <- c("z_AI","z_AI_pval")
							purity_ase <- tryCatch(cbind(df_info[c("RSID")], tst.bbreg.purity_res), error=function(e) data.frame(df_info[c("RSID")],z_AI=NA,z_AI_pval=NA))
							#setDF(purity_ase)
							purity_ase <- purity_ase %>% mutate_all(na_if,"")
							#tst.bbreg.purity_res = tryCatch(as.data.frame(tst.bbreg.purity)["cond",c(2,1)], error=function(e) data.frame(z_AI=NA,z_AI_pval=NA))
							#colnames(tst.bbreg.purity_res) <- c("z_AI","z_AI_pval")
							#purity_interaction <- tryCatch(cbind(df_info, tst.bbreg.purity_res , SAMPLE_NAME), error=function(e) data.frame(df_info,z_AI=NA,z_AI_pval=NA,SAMPLE_NAME))							
							#fwrite(purity_interaction, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".interactionASE_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".interactionASE_stratas_results.txt") )
							#message(paste0(s,": bbreg purity test finished"))
							rm(tst.bbreg.purity,tst.bbreg.purity_res)
						}
						if ( DO.CELLSPECIFIC ) {
								message("Doing Interaction ASE test")
									BBREG_list <- mclapply(cell_types, function(c){
										#c=cell_types[1]
										cf0 <- cell_pop_sub1[, c]
										all.cf0 = cf0[CUR.IND]
										all.cf0[is.na(all.cf0)] <-mean(all.cf0,na.rm=T)
										if (DO.RANKNORM) {
											all.cf0=rank.normalize(all.cf0)
										}
										tst.bbreg.cell = tryCatch(bbreg.cf0.test( CUR.REF, CUR.ALT, all.RHO, purity=all.NF_ase, cond=all.cf0), error=function(e) NA)
   										tst.bbreg.cell_res = tryCatch(as.data.frame(tst.bbreg.cell)["cond",c(2,1)], error=function(e) data.frame(z_AI=NA,z_AI_pval=NA))
   										colnames(tst.bbreg.cell_res) <- c("z_AI","z_AI_pval")
   										z_AI = tryCatch(data.frame(cell=c, z_AI=tst.bbreg.cell_res$z_AI,z_AI_pval=tst.bbreg.cell_res$z_AI_pval), error=function(e) data.frame(cell=c,z_AI=NA,z_AI_pval=NA))
										#tst.bbreg.cell = tryCatch(bbreg.cf0.test( CUR.REF, CUR.ALT, all.RHO, purity=(1-all.NF_ase), cond=all.cf0), error=function(e) NA)
   										#tst.bbreg.cell_res = tryCatch(as.data.frame(tst.bbreg.cell)["cond",c(2,1)], error=function(e) data.frame(z_AI=NA,z_AI_pval=NA))
   										#colnames(tst.bbreg.cell_res) <- c("z_AI","z_AI_pval")
   										#z_AI_TF = tryCatch(data.frame(z_AI_TF=tst.bbreg.cell_res$z_AI,z_AI_TF_pval=tst.bbreg.cell_res$z_AI_pval), error=function(e) data.frame(z_AI_TF=NA,z_AI_TF_pval=NA))
										#df <- tryCatch(cbind(df_info[c("RSID")], z_AI, z_AI_TF), error=function(e) data.frame(df_info[c("RSID")],cell=c,z_AI=NA,z_AI_pval=NA,z_AI_TF=NA,z_AI_TF_pval=NA))
										df <- tryCatch(cbind(df_info[c("RSID")], z_AI), error=function(e) data.frame(df_info[c("RSID")],cell=c,z_AI=NA,z_AI_pval=NA))
										#cat("inside BBREG loop",df$cell, df$z_AI, df$z_AI_pval, '\n' , sep='\t' )
										return(df)})
									z_BBREG_df <- ldply(BBREG_list, data.frame)
									rm(m1,m2,mat_s,GENO.H1_s,GENO.H2_s,BBREG_list,all.RHO,CUR.REF,CUR.ALT,CUR.IND,CUR.ID)
									#for QC while running uncomment
									#cat("BBREGtest",z_BBREG_df$z_AI[1], z_BBREG_df$z_AI_pval[1], '\n' , sep='\t' )
									#message(paste0(s,": bbreg test finished"))
   						}
   					}
				}
			}
			
			if (DO.CELLSPECIFIC & DO.BBREG ) {
				if (!exists('z_BBREG_df')) {
					#z_BBREG_df <- data.frame(df_info[c("RSID")],cell=NA,z_AI=NA,z_AI_pval=NA,z_AI_TF=NA,z_AI_TF_pval=NA)
					z_BBREG_df <- data.frame(df_info[c("RSID")],cell=NA,z_AI=NA,z_AI_pval=NA)
				}
				if (!exists('z_eqtl_df')) {
					#z_eqtl_df <- data.frame(df_info,cell=NA,z_eqtl=NA,z_eqtl_pval=NA,z_eqtl_TF=NA,z_eqtl_TF_pval=NA)
					z_eqtl_df <- data.frame(df_info,cell=NA,z_eqtl=NA,z_eqtl_pval=NA)
				}
				message(paste0(df_info[,"RSID"]," : Combining CF Interaction test"))
				sumQ <- merge(z_eqtl_df,z_BBREG_df, by=c("cell","RSID"),all=T)
			   	sumQ[mapply(is.infinite, sumQ)] <- NA
		    	#sumQ <- transform(sumQ, c=as.numeric(!is.na(z_eqtl)) + as.numeric(!is.na(z_AI)),c_TF=as.numeric(!is.na(z_eqtl_TF)) + as.numeric(!is.na(z_AI_TF)))
		    	sumQ <- transform(sumQ, c=as.numeric(!is.na(z_eqtl)) + as.numeric(!is.na(z_AI)))
		    	zs <- mclapply(split(sumQ,sumQ$cell), function(i){
		    		#i <- split(sumQ,sumQ$cell)[[4]]
				if(sum(i$c,na.rm=T)>0) {
		    		bothz <- c(i$z_eqtl,i$z_AI)
		    		stouffer <- transform(i, z_comb=sum(bothz,na.rm=T)/sqrt(c)) #3
		    		#stouffer <- transform(stouffer, z_comb_pval=2*pnorm(-abs(z_comb)))
		    		stouffer <- transform(stouffer, z_comb_pval=2*pnorm(-abs(z_comb)), sample=SAMPLE_NAME)
				} else {
				#stouffer <- transform(i, z_comb=NA, z_comb_pval=NA)
				stouffer <- transform(i, z_comb=NA, z_comb_pval=NA, sample=SAMPLE_NAME)
				}
				#if(sum(i$c_TF,na.rm=T)>0) {
				#bothz <- c(i$z_eqtl_TF,i$z_AI_TF)
		    		#stouffer <- transform(stouffer, z_comb_TF=sum(bothz,na.rm=T)/sqrt(c_TF)) #3
		    		#stouffer <- transform(stouffer, z_comb_TF_pval=2*pnorm(-abs(z_comb_TF)), sample=SAMPLE_NAME)
				#} else {
				#stouffer <- transform(stouffer, z_comb_TF=NA, z_comb_TF_pval=NA, sample=SAMPLE_NAME)
				#}   		
				return(stouffer)})
			sumQ <- ldply(zs, data.frame)
	    		#sumQ <- subset(sumQ[,-1], !is.na(z_comb) | !is.na(z_comb_TF))
	    		sumQ <- subset(sumQ[,-1], !is.na(z_comb) )
			if (dim(sumQ)[1] == 0) {
			#message("all na")
			} else {
			fwrite(sumQ, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".cfQTL.bbreg_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".cfQTL.bbreg_stratas_results.txt") )
			sumQ_header <- colnames(sumQ)
			}
			}
			if (DO.PURITY & DO.INTERACTION) {
				if (!exists('purity_ase')) {
					purity_ase <- data.frame(df_info[c("RSID")],z_AI=NA,z_AI_pval=NA)
				}
				if (!exists('purity_qtl')) {
					purity_qtl <- data.frame(df_info,z_eqtl=NA,z_eqtl_pval=NA)
				}
				message(paste0(df_info[,"RSID"]," : Combining NF Interaction test"))
				if ( "RSID" %in% colnames(purity_qtl) & "RSID" %in% colnames(purity_ase) ){
				sumQ <- merge(purity_qtl,purity_ase, by=c("RSID"),all=T)
				} else if ( "RSID" %in% colnames(purity_qtl) & !"RSID" %in% colnames(purity_ase) ){
				sumQ <- data.frame(purity_qtl,z_AI=NA,z_AI_pval=NA)
				} else if ( !"RSID" %in% colnames(purity_qtl) & "RSID" %in% colnames(purity_ase) ) {
				sumQ <- data.frame(df_info,z_eqtl=NA,z_eqtl_pval=NA,purity_ase[,"z_AI"],purity_ase[,"z_AI_pval"])
				} else {
				sumQ <- data.frame(df_info,z_eqtl=NA,z_eqtl_pval=NA,z_AI=NA,z_AI_pval=NA)
				}
		   	sumQ[mapply(is.infinite, sumQ)] <- NA
		    	sumQ <- transform(sumQ, c=as.numeric(!is.na(z_eqtl)) + as.numeric(!is.na(z_AI)))
			if(sum(sumQ$c,na.rm=T)>0) {
			bothz <- c(sumQ$z_eqtl,sumQ$z_AI)
		    	sumQ <- transform(sumQ, z_comb=sum(bothz,na.rm=T)/sqrt(c)) #3
		    	sumQ <- transform(sumQ, z_comb_pval=2*pnorm(-abs(z_comb)))
	    		sumQ <- subset(sumQ, !is.na(z_comb))
			fwrite(sumQ, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".purityQTL_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".purityQTL.stratas_results.txt") )
			} else {
				message("No NFqtl")
			}
			}
			if (DO.PRS & DO.BBREG) {
				if (!exists('PRS_ase')) {
					PRS_ase <- data.frame(df_info[c("RSID")],z_AI=NA,z_AI_pval=NA)
				}
				if (!exists('PRS_qtl')) {
					PRS_qtl <- data.frame(df_info,z_eqtl=NA,z_eqtl_pval=NA)
				}
				message(paste0(df_info[,"RSID"]," : Combining PRS Interaction test"))
				sumQ <- merge(PRS_qtl,PRS_ase, by=c("RSID"),all=T)
			
		   	sumQ[mapply(is.infinite, sumQ)] <- NA
		    	sumQ <- transform(sumQ, c=as.numeric(!is.na(z_eqtl)) + as.numeric(!is.na(z_AI)))
			if(sum(sumQ$c,na.rm=T)>0) {
			bothz <- c(sumQ$z_eqtl,sumQ$z_AI)
		    	sumQ <- transform(sumQ, z_comb=sum(bothz,na.rm=T)/sqrt(c)) #3
		    	sumQ <- transform(sumQ, z_comb_pval=2*pnorm(-abs(z_comb)))
	    		sumQ <- subset(sumQ, !is.na(z_comb))
			fwrite(sumQ, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".PRSQTL_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".PRSQTL.stratas_results.txt") )
			} else {
				message("No PRSqtl")
			}
			}
		})
	}
}

warnings()
