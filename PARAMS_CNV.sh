#!/bin/sh

# ---
# Parameters for this script

# Sample identifier from the ASVCF file
# Example: TCGA-CZ-5468
ID=$1

# CNV file (in TCGA/FireCloud format)
# Example file: KIRC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.seg.txt (unzip if gz)
# Example download: http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/KIRC/20160128/gdac.broadinstitute.org_KIRC.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.Level_3.2016012800.0.0.tar.gz
CNV=$2

# stratas directory
DIR=$3

# AS VCF prefix for each chromosome (assuming split by chromosome)
# {1-22}.vcf.gz will be appended to this to load the file
#VCF=$4

# ---

begin=$(date +"%s")
module load R/3.6.3

echo started counts
# Get the counts from each sample
# prev run bcftools query -f '%CHROM\t%END\t%ID\t%REF[\t%SAMPLE\t%GT\t%AS]\n' /data/gusev/USERS/ckal/BRCA/vcf/${CSV1}.reheadered.vcf | tr ',' '\t' > /data/gusev/USERS/ckal/BRCA/vcf/${CSV1}.counts

awk '{ print $1,$2,$6,$7,$8 }' ${DIR}/vcf/completed/${ID}.counts | awk 'BEGIN{ print "CHR POS HAP REF.READS ALT.READS" } $3 != "0|0" && $3 != "1|1" && ($4+$5) > 0 { print $0 }' | tr ' ' '\t' > ${DIR}/vcf/completed/${ID}.no0.counts

# Get the CNV segments
echo started CNV
cat ${DIR}/$CNV | grep $ID | awk 'BEGIN { print "ID CHR P0 P1 CNV" } { print $1,"chr"$2,$3,$4,$NF }' | tr ' ' '\t' > ${DIR}/$ID.cnv

# Compute CNVs for each region w/ more than 50 SNPs
echo started params
Rscript params.R --inp_counts ${DIR}/vcf/completed/${ID}.no0.counts --inp_cnv ${DIR}/$ID.cnv --out ${DIR}/$ID.par_region --min_cov 5 --id $ID --min_snps 50

# For any regions w/ fewer than 50 SNPs, rerun aggregating across regions
echo started aggregating
cat ${DIR}/$ID.par_region.local.params | awk 'NR == 1 || $6 == "NA"' | cut -f 1-5 > ${DIR}/$ID.rerun.cnv
echo ${DIR}/$ID.rerun.cnv
Rscript params.R --inp_counts ${DIR}/vcf/completed/$ID.no0.counts --inp_cnv ${DIR}/$ID.rerun.cnv --out ${DIR}/$ID.par_group --min_cov 5 --id $ID --group 10 --min_snps 50

# Merge the two outputs
echo started merging outputs
mv ${DIR}/$ID.par_region.global.params ${DIR}/$ID.global.params
cat ${DIR}/$ID.par_region.local.params | grep -v NA \
| cat - ${DIR}/$ID.par_group.local.params | sort -k2,2 -k3,3n | uniq > ${DIR}/$ID.local.params

# Clean up
echo started clean up
rm ${DIR}/$ID.par_region.* ${DIR}/$ID.par_group.* ${DIR}/$ID.rerun.cnv

termin=$(date +"%s")
difftimelps=$(($termin-$begin))
echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed for Script Execution. #Completed: ${set_num}"


# Some important variables to check (Can be removed later)
echo '---PROCESS RESOURCE LIMITS---'
ulimit -a
