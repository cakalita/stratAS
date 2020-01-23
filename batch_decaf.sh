#!/bin/bash
#SBATCH -n 10                # Number of cores (processes/threads?)
#SBATCH -N 1                # Ensure that all cores are on one machine; number of nodes
#SBATCH --mem=20000           # Memory pool for all cores (see also --mem-per-cpu)
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
echo started ${GENE_SET}
echo ${OUT_FOLDER}

#need matrix of variant info, genotypes, and ref/alt counts, all split by space and then I think it splits lines?
#this script runs stratas on a subset of data 

mkdir ${DATA}/stratas_results
mkdir ${DATA}/stratas_results/${OUT_FOLDER}

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
#also run for matrix file
#awk -v OFS='\t' '{print > ("KIRC.ALL.AS."$1".OUT.MAT.split")}' KIRC.ALL.AS.OUT.MAT.split

Rscript slurm_src/stratas_sasha_mc_gene_local.R \
	--input ${DATA}/stratas_input_files/${SAMPLE}.${CHR}.OUT.MAT.split \
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

#example run:
#for gene_set in ../../agusevlab/LNCaP/split_gene/*chr3*+([0-9]); do
#  sbatch --job-name=${gene_set} --export=BINOM="TRUE",INDIV="TRUE",LOCAL="KIRC.ALL.AS.CNVLOCAL",GLOBAL="KIRC.ALL.AS.nona.global.params",SAMPLE="KIRC.ALL.AS.01A",DATA="../../agusevlab/ckal/TCGA_vcf/KIRC",PEAK="LNCaP/split_gene/${gene_set}",GENE_SET="${gene_set}",CELL="../../agusevlab/ckal/AIsim/KIRC_TIMER.txt",OUT_FOLDER="test_batch/",EXPRESS="broad_rnaseq/KIRC.rnaseq.fixedsamples.full.txt" slurm_src/batch_decaf.sh 
#  sleep 1
#done
