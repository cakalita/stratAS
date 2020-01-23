#!/bin/sh
#SBATCH -n 10                # Number of cores (processes/threads?)
#SBATCH -N 1                # Ensure that all cores are on one machine; number of nodes
#SBATCH --mem=15000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o job_readouts/myoutput_%j.out  #"outFile"$1".txt"# File to which STDOUT will be written, %j inserts jobid
#SBATCH -e job_readouts/myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cakalita@gmail.com

begin=$(date +"%s")

#work from home folder
#stored big files in group
GROUP="../../agusevlab"

#load modules
module load R/3.5.2 
module load htslib

echo started ${GENE_SET}
echo ${OUT}

mkdir ${DATA}/stratas_results
mkdir ${DATA}/stratas_results/${OUT}

#NOTE that peak file needs to be preprocessed to have columns: CHR     P0      P1      NAME    CENTER

set_num=$(basename ${GENE_SET} | rev | cut -d"." -f1-2 | rev)
CHR=$(basename ${GENE_SET} | rev | cut -d"." -f2 | rev )

echo $set_num
echo $CHR

#run this before running script
#In agusevlab/LNCaP/split_gene
#for chr in {1..22}; do 
#to split gene file into chromosomes 
#	awk -v OFS='\t' '{print > ("../knownGene_withAlias.stratas."$1".txt")}' ../knownGene_withAlias.stratas.txt
#to split chr gene file into smaller parts
#	split -d -l 150 ../knownGene_withAlias.stratas.chr${chr}.txt knownGene_withAlias.stratas.chr${chr}.
#done
#for chr in {1..22}; do
#tabix -h ${DATA}/stratas_prep_files/${VCF}.vcf.gz chr${chr} | bgzip -c -f > ${DATA}/stratas_prep_files/${VCF}.chr${chr}.vcf.gz
#bcftools index ${DATA}/stratas_prep_files/${VCF}.chr${chr}.vcf.gz
#done

INP=${GENE_SET}
CHR=`cat $INP | tail -n+2 | head -n1 | awk '{ print $1 }'`
P0=`cat $INP | tail -n+2 | head -n1 | awk '{ print $2 - 100e3 }'`
P1=`cat $INP | tail -n1 | awk '{ print $3 + 100e3 }'`

OUT="$INP"

echo $OUT 
echo $CHR $P0 $P1

tabix ${DATA}/stratas_prep_files/${VCF}.${CHR}.vcf.gz "${CHR}:${P0}-${P1}" \
| cut -f 1-5,9- \
| tr ':' '\t' \
| awk '{ for(i=1;i<=5;i++) { printf "%s ",$i; } for(i=6;i<=10;i++) { if($i == "GT") gtnum=i-6; if($i=="AS") asnum=i-6; } for(i=11;i<=NF;i++) { if( (i-11) % 5 == asnum || (i-11) % 5 == gtnum ) printf " %s",$i; } print ""; }' \
| sed 's/[|,]/ /g' > $OUT.mat

if [ ! -s "$OUT.mat" ]; then
	echo $OUT.mat is empty
	exit
else

echo input is $OUT.mat
echo peaks is ${GENE_SET}
echo data path is ${DATA}/stratas_results/${OUT_FOLDER}

Rscript slurm_src/stratas_sasha_mc_gene_local.R \
	--input $OUT.mat \
	--global_param ${DATA}/stratas_input_files/${GLOBAL} \
	--local_param ${DATA}/CNV/${LOCAL} \
	--peaks ${GENE_SET} \
	--cell_specific TRUE \
	--cell_pop ${CELL} \
	--pheno FALSE \
	--gene_express ${DATA}/${EXPRESS} \
	--bbreg TRUE \
	--data_path ${DATA}/stratas_results/${OUT_FOLDER} \
	--betabinom ${BINOM} \
	--indiv ${INDIV} \
	--sample_name ${SAMPLE}.${set_num}
#--gene_name DPF3 #for running a single gene

termin=$(date +"%s")
difftimelps=$(($termin-$begin))
echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed for Script Execution. #Completed: ${set_num}"

fi

#example run: 
#for gene_set in ../../agusevlab/LNCaP/split_gene/*chr3*+([0-9]); do
#  sbatch --job-name=${gene_set} --export=BINOM="TRUE",INDIV="TRUE",LOCAL="KIRC.ALL.AS.CNVLOCAL",GLOBAL="KIRC.ALL.AS.nona.global.params",VCF="KIRC.ALL.AS.filtered_01A",SAMPLE="KIRC.ALL.AS.01A",DATA="../../agusevlab/ckal/TCGA_vcf/KIRC",PEAK="LNCaP/split_gene/${gene_set}",GENE_SET="${gene_set}",CELL="../../agusevlab/ckal/AIsim/KIRC_TIMER.txt",OUT_FOLDER="tabix_batch/",EXPRESS="broad_rnaseq/KIRC.rnaseq.fixedsamples.full.txt" slurm_src/tabix_batch_decaf.sh 
#  sleep 1
#done


