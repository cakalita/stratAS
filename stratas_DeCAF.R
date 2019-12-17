library(VGAM)
library(optparse)
library(plyr)
library(parallel)
library(dplyr)
library(data.table)

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
	make_option("--out", action="store", default=NA, type='character',
              help="Path to output [default: %default]"),
	make_option("--cell_specific", action="store_true", default=FALSE,
              help="Perform cell-type specific ASE test. [default: %default]"),
	make_option("--data_path", action="store", default=NA, type='character',
              help="Set output data path to store results. [default: %default]"),
	make_option("--pheno", action="store_true", default=TRUE,
              help="Preforms liklihood ratio test for chisq results between conditions. Curently binary only. [default: %default]"),
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

message("Files being read in")
peaks = read.table( opt$peaks , head=F , as.is=T, sep='\t')
#chr genes peaks file missing header
colnames(peaks) <- c("CHR","P0","P1","NAME","CENTER")
mat = read.table( opt$input , as.is=T )
cnv.all = read.table( opt$global_param , head=T ,as.is=T)
#cnv.all.both = read.table( opt$global_param , head=T ,as.is=T)
#cnv.all = cnv.all.both[grep("*[-]01A",cnv.all.both$ID),]
#the following is just until params
#samples file contains one sample ID per line 
#bcftools query -l KIRC.ALL.AS.merged.fixed.vcf.gz > samples_KIRC.ALL.AS.merged.fixed.full.txt
#next line subsets tumor samples from TCGA
#less samples_KIRC.ALL.AS.merged.fixed.full.txt | grep '.*01A$' > samples_KIRC.ALL.AS.merged.fixed.01A.txt

samples <- fread("../../agusevlab/ckal/TCGA_vcf/KIRC/stratas_prep_files/samples_KIRC.ALL.AS.merged.fixed.01A.txt", header=F)
cnv.all <- cnv.all[cnv.all$ID %in% samples$V1,]


if (DO.PHENO) {
	phe = read.table( opt$samples , head=T , as.is=T)
	message("Doing condition test")
}
if (DO.CELLSPECIFIC) {
	gene_express = read.table( opt$gene_express , head=T , as.is=T, sep='\t', check.names=FALSE) #tcga samples in header were getting converted from dast to dot
	cell_pop.both = read.table( opt$cell_pop , head=T , as.is=T)
	cell_pop = cell_pop.both[grep("*[-]01A",cell_pop.both$sample_ID),]
	#samples_counts = fread("stratas_prep_files/counts.samples.txt",header=F)
	message("Doing cell_specific test")
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

bbreg.test = function( ref , alt , rho , covar , cond=NULL ) {
    if ( !is.null(cond) ) {
    	if ( sum(!is.na(covar)) == 0 ) return( list( "pv" = c(NA,NA,NA) ) )
      df = data.frame( y=alt , n=ref+alt , cond = cond , covar = covar , rho=rho )      
    } else {
    	if ( sum(!is.na(covar)) == 0 ) return( list( "pv" = c(NA,NA) ) )
      df = data.frame( y=alt , n=ref+alt , covar = covar , rho=rho )
    }

    # test with a fixed overd parameter for each individual (for some reason this is very SLOW!)
    # reg = betabin( cbind( y , n - y  ) ~ 1 + cond , ~ phi.group , df , fixpar = list( 3:(n+2) , rho ) )
    df = df[ !is.na(df$covar) & df$n > 0, ]
    
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

    # need at least two samples to do regression with covariates

    if ( !is.null(cond) ) {
      # need at least two samples in both conditions
	statement= nrow(df) < 2 || cor(df$covar,df$cond,use = "complete.obs") == 1 || min(rle(sort(df$cond))$len) < 2 || sd(df$covar,na.rm=T) == 0 || sd(df$cond,na.rm=T) == 0
      if ( statement ) return( list( "pv" = c(NA,NA,NA) ) )
      else reg = betabin( cbind( y , n - y  ) ~ 1 + cond + covar , ~ phi.q , df , fixpar = list( 4:(nq+3) , rho.q ) )
    } else {
      if ( nrow(df) < 2 || sd(df$covar,na.rm=T) == 0 ) return( list( "pv" = c(NA,NA) ) )
      else reg = betabin( cbind( y , n - y  ) ~ 1 + covar , ~ phi.q , df , fixpar = list( 3:(nq+2) , rho.q ) )
    }
    
    zscores = coef(reg) / sqrt(diag(vcov( reg )))
    pvals = 2*(pnorm( abs(zscores) , lower.tail=F))
    opt = list( "pv" = pvals, "z" = zscores )
  return( opt )
}


bbreg.cf0.test = function( ref , alt , rho , covar , cond=NULL ) {
	sum.covar = sum(!is.na(covar))
    if ( !is.null(cond) ) {
    	if ( sum.covar == 0 ) return( list( "pv" = c(NA,NA,NA) ) )
      df = data.frame( y=alt , n=ref+alt , cond = cond , covar = covar , rho=rho )      
    } else {
    	if ( sum.covar == 0 ) return( list( "pv" = c(NA,NA) ) )
      df = data.frame( y=alt , n=ref+alt , covar = covar , rho=rho )
    }

    # test with a fixed overd parameter for each individual (for some reason this is very SLOW!)
    # reg = betabin( cbind( y , n - y  ) ~ 1 + cond , ~ phi.group , df , fixpar = list( 3:(n+2) , rho ) )
    df = df[ !is.na(df$covar) & df$n > 0, ]
    
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

    # need at least two samples to do regression with covariates

    if ( !is.null(cond) ) {
      # need at least two samples in both conditions
	length_condition = min(rle(sort(df$cond))$len) < 2
	statement= nrow(df) < 2 || cor(df$covar,df$cond,use = "complete.obs") == 1 || sd(df$covar,na.rm=T) == 0 || sd(df$cond,na.rm=T) == 0
      if ( statement) return( list( "pv" = c(NA,NA,NA) ) )
      else reg = betabin( cbind( y , n - y  ) ~ 1 + cond + covar , ~ phi.q , df , fixpar = list( 4:(nq+3) , rho.q ) )
    } else {
      if ( nrow(df) < 2 || sd(df$covar,na.rm=T) == 0 ) return( list( "pv" = c(NA,NA) ) )
      else reg = betabin( cbind( y , n - y  ) ~ 1 + covar , ~ phi.q , df , fixpar = list( 3:(nq+2) , rho.q ) )
    }
    
    zscores = coef(reg) / sqrt(diag(vcov( reg )))
    pvals = 2*(pnorm( abs(zscores) , lower.tail=F))
    opt = list( "pv" = pvals, "z" = zscores )
  return( opt )
}

cur.chr = unique(mat[,1])
m = match(peaks$CHR , cur.chr )
peaks = peaks[!is.na(m),]

if ( !is.na(opt$local_param) )  {
	cnv.local = read.table( opt$local_param , head=T ,as.is=T)
	#cnv.local.both = read.table( opt$local_param , head=T ,as.is=T)
	#cnv.local = cnv.local.both[grep("*[-]01A",cnv.local.both$ID),]
	m = match(cnv.local$CHR , cur.chr )
	cnv.local = cnv.local[!is.na(m),]
}

N = (ncol(mat) - 5)/4
M = nrow(mat)

HAPS = list()
GENO.H1 = matrix( 0 , nrow=M , ncol=N )
GENO.H2 = matrix( 0 , nrow=M , ncol=N )

if (DO.PHENO) {
	PHENO = phe$CONDITION
	m = match( phe$ID , cnv.all$ID )
	cnv.all = cnv.all[m,]
	RHO.ALL = cnv.all$PHI
} else {
	RHO.ALL = cnv.all$PHI
}

for ( h in 1:2 ) {
	HAPS[[h]] = matrix( 0 , nrow=M , ncol=N )
}

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

options( digits = 4 )

RHO = RHO.ALL
RHO[ RHO > MAX.RHO ] = NA

if (DO.CELLSPECIFIC) {
	gene = merge(peaks, gene_express, by="NAME")
}

COL.HEADER = c("CHR","POS","RSID","P0","P1","NAME","CENTER","N.HET","N.READS")
COL.HEADER.PHENO = c("ALL.AF","ALL.BBINOM.P","C0.AF","C0.BBINOM.P","C1.AF","C1.BBINOM.P","DIFF.BBINOM.P")
COL.HEADER.BINOM = c("ALL.BINOM.P","ALL.C0.BINOM.P","ALL.C1.BINOM.P","FISHER.OR","FISHER.DIFF.P")
COL.HEADER.PHENO.INDIV = c("IND.C0","IND.C0.COUNT.REF","IND.C0.COUNT.ALT","IND.C1","IND.C1.COUNT.REF","IND.C1.COUNT.ALT")
COL.HEADER.INDIV = c("IND","IND.COUNT.REF","IND.COUNT.ALT")
COL.HEADER.BBREG = c("C0.BBREG.P","C0.CNV.BBREG.P","C1.BBREG.P","C1.CNV.BBREG.P","ALL.BBREG.P","DIFF.BBREG.P","DIFF.CNV.BBREG.P")
COL.HEADER.BETABINOM = c("ALL.AF","ALL.BBINOM.P")

HEAD = COL.HEADER
#if ( DO.PHENO ) HEAD = c(COL.HEADER,COL.HEADER.PHENO)
#if ( DO.BINOM ) HEAD = c(COL.HEADER,COL.HEADER.BINOM)
#if ( DO.BBREG ) HEAD = c(COL.HEADER,COL.HEADER.BBREG)
#if ( DO.INDIV ) HEAD = c(COL.HEADER,COL.HEADER.INDIV)
#if ( DO.BETABINOM ) HEAD = c(COL.HEADER,COL.HEADER.BETABINOM)
#cat( HEAD , sep='\t')
#cat('\n')

#for testing only
if ( !is.na(opt$gene_name) ) {
peaks = peaks[peaks$NAME==GENE_NAME,]
gene = gene[gene$NAME==GENE_NAME,]
}
for ( p in 1:nrow(peaks) ) {
#peak_list <- mclapply(1:nrow(peaks), function(p) {
	#p=1
	#browser()
	if (DO.CELLSPECIFIC) {
		Y.total = gene[p,-c(1:5)]
	}
	
	#cur_peak <- subset(mat, V1==peaks$CHR[p])
	#cur_peak <- dplyr::filter(cur_peak, between(cur_peak$V2, peaks$P0[p], peaks$P1[p]))
	#cur <- ifelse(mat$V3 %in% cur_peak$V3, TRUE, FALSE)
	cur = mat[,1] == peaks$CHR[p] & mat[,2] >= peaks$P0[p] & mat[,2] <= peaks$P1[p]

	if(sum(cur) == 0 ) next()

	if ( !is.na(opt$local_param) )  {
		cnv.local.cur = cnv.local[cnv.local$ID %in% cnv.all$ID,]
		#cnv.local.cur = cnv.local[cnv.local$ID %in% samples_counts$V1,]
		cur.cnv = cnv.local.cur[ cnv.local.cur$CHR == peaks$CHR[p] & cnv.local.cur$P0 < peaks$P0[p] & cnv.local.cur$P1 > peaks$P1[p] , ]
		#cur.cnv = cnv.local[ cnv.local$CHR == peaks$CHR[p] & cnv.local$P0 < peaks$P0[p] & cnv.local$P1 > peaks$P1[p] , ]
	
	if (DO.PHENO) {
			m = match( phe$ID , cur.cnv$ID )
			cur.cnv = cur.cnv[m,]
		}
		cnv.a <- merge(cnv.all,cur.cnv, by="ID", all=T)
		cnv.a <- subset(cnv.a, ID %in% samples$V1)
		#remove the following line when going to full individual counts
		cnv.a <- transform(cnv.a, PHI=ifelse(is.na(PHI.y), PHI.x, PHI.y))
		#cnv.a <- transform(cur.cnv, PHI=ifelse(is.na(PHI), 0.0263, PHI))
		RHO = cnv.a$PHI.y
		#RHO = cur.cnv$PHI
		RHO[ RHO > MAX.RHO ] = NA
		COVAR = cnv.a$CNV
	}

	# collapse reads at this peak
	cur.h1 = vector()
	cur.h2 = vector()
	cur.i = vector()
	
	for ( i in 1:N ) {
		reads.keep = GENO.H1[cur,i] != GENO.H2[cur,i] & HAPS[[1]][cur,i] >= MIN.COV & HAPS[[2]][cur,i] >= MIN.COV
		reads.keep[reads.keep] <- !c(FALSE, diff(mat[cur, 2][reads.keep]) < opt$exclude)
		cur.h1 = c( cur.h1 , (HAPS[[1]][cur,i])[reads.keep] )
		cur.h2 = c( cur.h2 , (HAPS[[2]][cur,i])[reads.keep] )
		cur.i = c( cur.i , rep( i , sum(reads.keep)) )
	}
	
	if ( length(unique(cur.i)) > MIN.MAF*N && sum(cur.h1) + sum(cur.h2) > 0 ) {
		# test all nearby SNPs

		if( PAR.WIN == -1 ) {
			cur.snp = mat[,1] == peaks$CHR[p] & mat[,2] >= peaks$P0[p] & mat[,2] <= peaks$P1[p]
		} else {
			cur.snp = mat[,1] == peaks$CHR[p] & mat[,2] >= peaks$CENTER[p] - PAR.WIN & mat[,2] <= peaks$CENTER[p] + PAR.WIN
		}

		#for ( s in which(cur.snp) ) {
		snp_list <- lapply(which(cur.snp), function(s){
		#s = 1
		#s = 117694
			# restrict to hets
			HET = GENO.H1[s,] != GENO.H2[s,] & !is.na(RHO)
	if ( !is.na(opt$local_param) )  {
			IDS = cnv.a$ID
			#IDS = cur.cnv$ID
			#IDS = samples_counts$V1
	} else {
			IDS = cnv.all$ID
	}
			# collect REF/ALT heterozygous haplotypes
				m1 = match( cur.i , which(HET & GENO.H1[s,] == 0) )
				m2 = match( cur.i , which(HET & GENO.H2[s,] == 0) )
				CUR.REF = c( cur.h1[ !is.na(m1) ] , cur.h2[ !is.na(m2) ])
				CUR.ALT = c( cur.h2[ !is.na(m1) ] , cur.h1[ !is.na(m2) ])
				CUR.IND = c( cur.i[ !is.na(m1) ] , cur.i[ !is.na(m2)] )
				CUR.ID = unique(IDS[CUR.IND]) #this is on the theory that CUR.IND is vector of Nth individuals representing that snp
				df_ind <- cbind(paste(CUR.IND,collapse=',') , paste(CUR.REF,collapse=',') , paste(CUR.ALT,collapse=',') , paste(CUR.ID,collapse=','))
				df_info <- cbind(mat[s,1:3] , peaks[p,c("P0","P1","NAME","CENTER")] , sum(HET) , sum(CUR.REF) + sum(CUR.ALT))
		          	colnames(df_info) <- COL.HEADER				
				if (DO.CELLSPECIFIC) {
				SNP_geno=GENO.H1[s,] + GENO.H2[s,]
				#cell_pop_sub1=subset(cell_pop, sample_ID %in% IDS)
				cell_pop_nodup <- cell_pop[!duplicated(cell_pop$sample_ID),]
				cnv.ids <- data.frame(ID=cnv.all[,"ID"])
				cell_pop_sub1 <- merge(cnv.ids, cell_pop_nodup, by.x="ID", by.y="sample_ID",all.x=T)
				eqtl_IDset <- ifelse(cell_pop_sub1$ID %in% colnames(Y.total), TRUE, FALSE)
				cell_pop_sub <- subset(cell_pop_sub1, ID %in% colnames(Y.total))
				cur.y.total = colnames(Y.total) %in% IDS
				Y.total.ind <- unname(unlist(Y.total[cur.y.total]))
				cell_types <- colnames(cell_pop_sub[,-1])
				SNP <- SNP_geno[eqtl_IDset]
 				COVAR_sub <- COVAR[eqtl_IDset]
				lm_eqtl_vanilla = tryCatch(lm(Y.total.ind ~ SNP + COVAR_sub), error=function(e) NA)
       				lm_eqtl_vanilla_df = as.data.frame(summary(lm_eqtl_vanilla)$coefficients)
   				lm_eqtl_vanilla_interaction = tryCatch(lm_eqtl_vanilla_df["SNP",c(3,4)], error=function(e) data.frame(z_eqtl_vanilla=NA,eqtl_pval_vanilla=NA))
    				colnames(lm_eqtl_vanilla_interaction) <- c("z_eqtl_vanilla","eqtl_pval_vanilla")
    				z_eqtl_vanilla = data.frame(z_eqtl_vanilla=lm_eqtl_vanilla_interaction$z_eqtl_vanilla,eqtl_pval_vanilla=2*pnorm(-abs(lm_eqtl_vanilla_interaction$z_eqtl_vanilla)))
				write.table(cbind(df_info,z_eqtl_vanilla, SAMPLE_NAME), quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".eqtl_vanilla_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".eqtl_vanilla_stratas_results.txt") )
				eqtl_list <- lapply(cell_types, function(c){
					#c=cell_types[1]
       				#cf0_L <- cell_pop_sub[, ..c]
       				#cf0 <- unname(unlist(cf0_L))
       				cf0 <- cell_pop_sub[, c]
 				lm_eqtl = tryCatch(lm(Y.total.ind ~ cf0*SNP + SNP + COVAR_sub), error=function(e) NA)
       				lm_eqtl_df = as.data.frame(summary(lm_eqtl)$coefficients)
   				lm_interaction = tryCatch(lm_eqtl_df["cf0:SNP",c(3,4)], error=function(e) data.frame(z_eqtl=NA,pval=NA))
    				colnames(lm_interaction) <- c("z_eqtl","pval")
    				z_eqtl = data.frame(cell=c, z_eqtl=lm_interaction$z_eqtl,z_eqtl_pval=2*pnorm(-abs(lm_interaction$z_eqtl)))
				df <- cbind(df_info, z_eqtl)
    				return(df)})
				z_eqtl_df <- ldply(eqtl_list, data.frame)
 				}	
			if ( sum(HET) > MIN.HET*N ) {

				if ( sum(CUR.REF) + sum(CUR.ALT) > 0 ) {
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

						if (DO.PHENO) {
						m = !is.na(match( CUR.IND , which(CUR.PHENO==0) ))
						CUR.REF.C0 = CUR.REF[m]
						CUR.ALT.C0 = CUR.ALT[m]
						CUR.IND.C0 = CUR.IND[m]

						m = !is.na(match( CUR.IND , which(CUR.PHENO==1) ))
						CUR.REF.C1 = CUR.REF[m]
						CUR.ALT.C1 = CUR.ALT[m]
						CUR.IND.C1 = CUR.IND[m]		
						
						df_ind <- cbind(paste(CUR.IND.C0,collapse=',') , paste(CUR.REF.C0,collapse=',') , paste(CUR.ALT.C0,collapse=',') , paste(CUR.IND.C1,collapse=',') , paste(CUR.REF.C1,collapse=',') , paste(CUR.ALT.C1,collapse=',') , paste(CUR.ID,collapse=','))

						# --- perform beta-binomial test
						tst.bbinom.C0 = bbinom.test( CUR.REF.C0 , CUR.ALT.C0 , RHO[CUR.IND.C0] )
						tst.bbinom.C1 = bbinom.test( CUR.REF.C1 , CUR.ALT.C1 , RHO[CUR.IND.C1] )
						tst.bbinom.ALL = bbinom.test( CUR.REF , CUR.ALT , RHO[CUR.IND] )
						lrt.BOTH = 2 * (tst.bbinom.C0$objective + tst.bbinom.C1$objective - tst.bbinom.ALL$objective)
						pv.BOTH = pchisq( abs(lrt.BOTH) , df=1 , lower.tail=F )			

						# --- print main output
						#cat( unlist(mat[s,1:3]) , unlist(peaks[p,c("P0","P1","NAME","CENTER")]) , sum(HET) , sum(CUR.REF) + sum(CUR.ALT) , tst.bbinom.ALL$min , tst.bbinom.ALL$pv , tst.bbinom.C0$min , tst.bbinom.C0$pv , tst.bbinom.C1$min , tst.bbinom.C1$pv , pv.BOTH , SAMPLE_NAME , '\n' , sep='\t', file=paste0(DATA_PATH,SAMPLE_NAME,".bbinom_stratas_results.txt"), append=TRUE )
						HEAD = c(COL.HEADER,COL.HEADER.PHENO)
						df <- cbind(df_info, tst.bbinom.ALL$min , tst.bbinom.ALL$pv , tst.bbinom.C0$min , tst.bbinom.C0$pv , tst.bbinom.C1$min , tst.bbinom.C1$pv , pv.BOTH , SAMPLE_NAME)
						colnames(df) <- HEAD
						write.table(df, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".bbinomLRT_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".bbinomLRT_stratas_results.txt") )
							# --- print individual counts
								if ( DO.INDIV ) {
									HEAD = c(HEAD,COL.HEADER.PHENO.INDIV)
								  	#cat( "" , paste(CUR.IND.C0,collapse=',') , paste(CUR.REF.C0,collapse=',') , paste(CUR.ALT.C0,collapse=',') , paste(CUR.IND.C1,collapse=',') , paste(CUR.REF.C1,collapse=',') , paste(CUR.ALT.C1,collapse=',') , paste(CUR.ID,collapse=','), SAMPLE_NAME , '\n' , sep='\t', file=paste0(DATA_PATH,SAMPLE_NAME,".bbinom.ind._stratas_results.txt"), append=TRUE )
									df <- cbind(df, df_ind)
									colnames(df) <- HEAD
									write.table(df, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".bbinomLRT.ind_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".bbinomLRT.ind_stratas_results.txt") )
								}

							# --- perform binomial test
							if ( DO.BINOM ) {
								HEAD = c(COL.HEADER,COL.HEADER.BINOM)
								tst.binom = binom.test( sum(CUR.ALT),sum(CUR.REF)+sum(CUR.ALT) )
								if ( sum(CUR.REF.C0)+sum(CUR.ALT.C0) > 0 ) tst.binom.C0 = binom.test( sum(CUR.ALT.C0),sum(CUR.REF.C0)+sum(CUR.ALT.C0) ) else tst.binom.C0 = list("p.value"=NA)
								if ( sum(CUR.REF.C1)+sum(CUR.ALT.C1) > 0 ) tst.binom.C1 = binom.test( sum(CUR.ALT.C1),sum(CUR.REF.C1)+sum(CUR.ALT.C1) ) else tst.binom.C1 = list("p.value"=NA)
								# --- perform fisher's exact test between conditions					
								tst.fisher = fisher.test( cbind( c(sum(CUR.REF.C0),sum(CUR.ALT.C0)) , c(sum(CUR.REF.C1),sum(CUR.ALT.C1)) ) )
								#cat( "" , tst.binom$p.value ,  tst.binom.C0$p.value , tst.binom.C1$p.value , tst.fisher$est , tst.fisher$p.value , SAMPLE_NAME , '\n' , sep='\t', file=paste0(DATA_PATH,SAMPLE_NAME,".binom_stratas_results.txt"), append=TRUE )
								df <- cbind(df_info, tst.binom$p.value ,  tst.binom.C0$p.value , tst.binom.C1$p.value , tst.fisher$est , tst.fisher$p.value , SAMPLE_NAME)
								colnames(df) <- HEAD
								write.table(df, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".binom_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".binom_stratas_results.txt") )
														# --- print individual counts
								if ( DO.INDIV ) {
									HEAD = c(HEAD,COL.HEADER.PHENO.INDIV)
								  	#cat( "" , paste(CUR.IND.C0,collapse=',') , paste(CUR.REF.C0,collapse=',') , paste(CUR.ALT.C0,collapse=',') , paste(CUR.IND.C1,collapse=',') , paste(CUR.REF.C1,collapse=',') , paste(CUR.ALT.C1,collapse=',') , paste(CUR.ID,collapse=','), SAMPLE_NAME , '\n' , sep='\t', file=paste0(DATA_PATH,SAMPLE_NAME,".bbinom.ind._stratas_results.txt"), append=TRUE )
									df <- cbind(df, df_ind)
									colnames(df) <- HEAD
									write.table(df, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".binom.ind_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".binom.ind_stratas_results.txt") )
								}
							}
							
							# --- perform beta-binom regression
							if ( DO.BBREG ) {
								HEAD = c(COL.HEADER,COL.HEADER.BBREG)
							  	tst.bbreg.c0 = bbreg.test( CUR.REF.C0 , CUR.ALT.C0 , RHO[ CUR.IND.C0 ] , COVAR[CUR.IND.C0] )
							  	tst.bbreg.c1 = bbreg.test( CUR.REF.C1 , CUR.ALT.C1 , RHO[ CUR.IND.C1 ] , COVAR[CUR.IND.C1] )
							  	tst.bbreg = bbreg.test( CUR.REF , CUR.ALT , RHO[ CUR.IND ] , COVAR[CUR.IND] , CUR.PHENO[CUR.IND] )
							  
							  
							  #cat( "" , tst.bbreg.c0$pv , tst.bbreg.c1$pv , tst.bbreg$pv , SAMPLE_NAME , '\n' , sep='\t', file=paste0(DATA_PATH,SAMPLE_NAME,"_stratas_results.txt"), append=TRUE )
								df <- cbind(df_info, tst.bbreg.c0$pv , tst.bbreg.c1$pv , tst.bbreg$pv , SAMPLE_NAME)
								colnames(df) <- HEAD
								write.table(df, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".bbreg_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".bbreg_stratas_results.txt") )
													# --- print individual counts
								if ( DO.INDIV ) {
									HEAD = c(HEAD,COL.HEADER.PHENO.INDIV)
								  	#cat( "" , paste(CUR.IND.C0,collapse=',') , paste(CUR.REF.C0,collapse=',') , paste(CUR.ALT.C0,collapse=',') , paste(CUR.IND.C1,collapse=',') , paste(CUR.REF.C1,collapse=',') , paste(CUR.ALT.C1,collapse=',') , paste(CUR.ID,collapse=','), SAMPLE_NAME , '\n' , sep='\t', file=paste0(DATA_PATH,SAMPLE_NAME,".bbinom.ind._stratas_results.txt"), append=TRUE )
									df <- cbind(df, df_ind)
									colnames(df) <- HEAD
									write.table(df, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".bbreg.ind_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".bbreg.ind_stratas_results.txt") )
								}
							}
							#cat('\n')
						}
						if ( DO.BETABINOM ) {
						# --- preform beta-binomial test without LRT to compare conditions
							HEAD = c(COL.HEADER,COL.HEADER.BETABINOM,"SAMPLE_NAME")
							tst.bbinom.ALL = tryCatch(bbinom.test( CUR.REF , CUR.ALT , RHO[CUR.IND] ), error=function(e) NA)
							#cat( unlist(mat[s,1:3]) , unlist(peaks[p,c("P0","P1","NAME","CENTER")]) , sum(HET) , sum(CUR.REF) + sum(CUR.ALT) , tst.bbinom.ALL$min , tst.bbinom.ALL$pv , '\n' , sep='\t' )
							#cat( unlist(mat[s,1:3]) , unlist(peaks[p,c("P0","P1","NAME","CENTER")]) , sum(HET) , sum(CUR.REF) + sum(CUR.ALT) , tst.bbinom.ALL$min , tst.bbinom.ALL$pv , SAMPLE_NAME , '\n' , sep='\t', file=paste0(DATA_PATH,SAMPLE_NAME,"_stratas_results.txt"), append=TRUE )
							df <- tryCatch(cbind(df_info , tst.bbinom.ALL$min , tst.bbinom.ALL$pv , SAMPLE_NAME), error=function(e) data.frame(df_info,AF=tst.bbinom.ALL$min,pval=tst.bbinom.ALL$pv,SAMPLE_NAME))
							colnames(df) <- HEAD
							write.table(df, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".bbinom_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".bbinom_stratas_results.txt") )
						# --- print individual counts
							if ( DO.INDIV ) {
								HEAD = c(HEAD,COL.HEADER.INDIV)
							  	#cat( "" , paste(CUR.IND.C0,collapse=',') , paste(CUR.REF.C0,collapse=',') , paste(CUR.ALT.C0,collapse=',') , paste(CUR.IND.C1,collapse=',') , paste(CUR.REF.C1,collapse=',') , paste(CUR.ALT.C1,collapse=',') , paste(CUR.ID,collapse=','), SAMPLE_NAME , '\n' , sep='\t', file=paste0(DATA_PATH,SAMPLE_NAME,".bbinom.ind._stratas_results.txt"), append=TRUE )
								df <- cbind(df, df_ind)
								colnames(df) <- HEAD
								write.table(df, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".bbinom.ind_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".bbinom.ind_stratas_results.txt") )
							}
						}
						if ( DO.BBREG & !DO.PHENO) {
							all.COVAR = COVAR[CUR.IND] 
							#all.COVAR[is.na(all.COVAR)] <- 0
							tst.bbreg.vanilla = tryCatch(bbreg.test( CUR.REF , CUR.ALT , RHO[ CUR.IND ] , all.COVAR ), error=function(e) NA)
							#cat( "" , tst.bbreg.c0$pv , tst.bbreg.c1$pv , tst.bbreg$pv , SAMPLE_NAME , '\n' , sep='\t', file=paste0(DATA_PATH,SAMPLE_NAME,"_stratas_results.txt"), append=TRUE )
							tst.bbreg.vanilla_res = tryCatch(as.data.frame(tst.bbreg.vanilla)["(Intercept)",c(2,1)], error=function(e) data.frame(z_AI=NA,z_AI_pval=NA))
							df <- tryCatch(cbind(df_info, tst.bbreg.vanilla_res , SAMPLE_NAME), error=function(e) data.frame(df_info,z_AI=NA,z_AI_pval=NA,SAMPLE_NAME))
							write.table(df, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".bbreg_vanilla_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".bbreg_vanilla_stratas_results.txt") )
						}
						if (DO.CELLSPECIFIC) {
							all.COVAR = COVAR[CUR.IND] 
							#all.COVAR[is.na(all.COVAR)] <- 0
							all.RHO = RHO[CUR.IND]
							#AI_list <- mclapply(cell_types, function(c){
						#		#c=cell_types[1]
	       					#		cf0_L <- cell_pop_sub1[, ..c]
						#		cf0 <- unname(unlist(cf0_L))
						#		all.cf0 = cf0[CUR.IND]
    						#	reg <- tryCatch(VGAM::vglm(cbind(CUR.REF,CUR.ALT) ~ all.cf0, betabinomial, irho=mean(all.RHO,na.rm=T), trace = FALSE), error=function(e) NA)
   						#		reg_res = tryCatch(as.data.frame(t(coef(summaryvglm(reg))["all.cf0",c(3,4)])), error=function(e) data.frame(z_AI=NA,z_AI_pval=NA))
   						#		colnames(reg_res) <- c("z_AI","z_AI_pval")
   						#		z_AI = data.frame(cell=c, z_AI=reg_res$z_AI,z_AI_pval=reg_res$z_AI_pval)
						#		df <- cbind(df_info[c("RSID")], z_AI)
    						#return(df)})
						#	z_AI_df <- ldply(AI_list, data.frame)

								if ( DO.BBREG ) {
									BBREG_list <- mclapply(cell_types, function(c){
										#c=cell_types[1]
										cf0_L <- cell_pop_sub1[, c]
										rownames(cf0_L) <- cell_pop_sub1$V1
										cell_pop_sub2<- na.omit(as.data.frame(cf0_L))
										cf0 <-cell_pop_sub2[,1] 
										all.cf0 = cf0[CUR.IND]
										all.cf0[is.na(all.cf0)] <-0
										#all.RHO[is.na(all.RHO)] <- 0.0263
										tst.bbreg.cell = tryCatch(bbreg.cf0.test( CUR.REF, CUR.ALT, all.RHO, covar=all.COVAR, cond=all.cf0), error=function(e) NA)
   										tst.bbreg.cell_res = tryCatch(as.data.frame(tst.bbreg.cell)["cond",c(2,1)], error=function(e) data.frame(z_AI=NA,z_AI_pval=NA))
   										colnames(tst.bbreg.cell_res) <- c("z_AI","z_AI_pval")
   										z_AI = data.frame(cell=c, z_AI=tst.bbreg.cell_res$z_AI,z_AI_pval=tst.bbreg.cell_res$z_AI_pval)
										df <- cbind(df_info[c("RSID")], z_AI)
    								return(df)})
									z_BBREG_df <- ldply(BBREG_list, data.frame)
									#for QC while running uncomment
									#cat("BBREGtest",z_BBREG_df$z_AI, z_BBREG_df$z_AI_pval, '\n' , sep='\t' )
   								}
   						}
   					}
				}
			}
			if (DO.CELLSPECIFIC) {
			if (DO.BBREG) {
			#z_BBREG_df <- tryCatch(data.frame(z_BBREG_df), error=function(e) data.frame(RSID=paste(df_info$RSID), z_AI=NA,z_AI_pval=NA,cell=NA))
			if (exists('z_BBREG_df') & exists('z_eqtl_df')){
			sumQ <- merge(z_eqtl_df,z_BBREG_df, by=c("cell","RSID"),all=T)
			} else if (exists('z_eqtl_df') & !exists('z_BBREG_df')) {
			sumQ <- data.frame(z_eqtl_df[,c(3,10,1:2,4:9,11:12)], z_AI=NA,z_AI_pval=NA)
			} else if (exists('z_BBREG_df') & !exists('z_eqtl_df')) {
			sumQ <- data.frame(cell=z_BBREG_df$cell, df_info[,c(3,1:2,4:9)], z_BBREG_df[,-2])
			} else {
				sumQ <-data.frame(cell=NA,rsID=NA,df_info[,-3],z_eqtl=NA,z_eqtl_pval=NA,z_AI=NA,z_AI_pval=NA )
				message("No ASE")
			}
		    sumQ[mapply(is.infinite, sumQ)] <- NA
		    sumQ <- transform(sumQ, c=ifelse(is.na(z_eqtl) & is.na(z_AI), 0,
		        ifelse(is.na(z_eqtl) & !is.na(z_AI),1,
		            ifelse(!is.na(z_eqtl) & is.na(z_AI),1, 2))))
		    zs <- mclapply(split(sumQ,sumQ$cell), function(i){
		    	#i <- split(sumQ,sumQ$cell)[[1]]
		    	zs <- c(i$z_eqtl,i$z_AI)
		    	sumQ <- transform(i, z_comb=sum(zs,na.rm=T)/sqrt(c)) #3
		    	sumQ <- transform(sumQ, z_comb_pval=2*pnorm(-abs(z_comb)), sample=SAMPLE_NAME)
		    	return(sumQ)})
			sumQ <- ldply(zs, data.frame)
	    		sumQ <- subset(sumQ[,-1], !is.na(z_comb))
			#for QC while running uncomment
			#cat("sumQ",sumQ$z_AI, sumQ$z_AI_pval, '\n' , sep='\t' )
			write.table(sumQ, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".cfQTL.bbreg_stratas_results_new.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".cfQTL.bbreg_stratas_results.txt") )
												# --- print individual counts
				if ( DO.INDIV ) {
				colnames(df_ind) <- COL.HEADER.INDIV
				df <- cbind(sumQ, df_ind)
				write.table(df, quote=F , row.names=F , sep='\t' , col.names=!file.exists(paste0(DATA_PATH,SAMPLE_NAME,".cfQTL.bbreg.ind_stratas_results.txt")) , append=TRUE, file=paste0(DATA_PATH,SAMPLE_NAME,".cfQTL.bbreg.ind_stratas_results.txt") )
				}
			}
			}
		if (DO.CELLSPECIFIC & DO.BBREG)	return(sumQ)
		})
	}
}
#if (DO.CELLSPECIFIC & DO.BBREG) {
#snps_df <- ldply(snp_list, data.frame)
#write.table(snps_df, quote=F , row.names=F , sep='\t' , file=paste0(DATA_PATH,SAMPLE_NAME,".cfQTL.bbreg_stratas_results_atend_new.txt") )
#}



