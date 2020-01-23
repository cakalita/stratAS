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
echo started ${CHR}

#need matrix of variant info, genotypes, and ref/alt counts, all split by space and then I think it splits lines?

mkdir ${DATA}/stratas_results
mkdir ${DATA}/stratas_results/${OUT_FOLDER}

#NOTE that peak file needs to be preprocessed to have columns: CHR     P0      P1      NAME    CENTER

echo ${DATA}/stratas_input_files/${SAMPLE}.chr${CHR}.OUT.MAT.split
echo ${DATA}/stratas_input_files/${GLOBAL}
echo ${DATA}/CNV/${LOCAL}
echo ${GROUP}/LNCaP/knownGene_withAlias.stratas.chr${CHR}.txt 
echo ${CELL}
echo ${DATA}/${EXPRESS}
echo ${DATA}/stratas_results/${OUT_FOLDER}

Rscript slurm_src/stratas_DeCAF.R \
	--input ${DATA}/stratas_input_files/${SAMPLE}.chr${CHR}.OUT.MAT.split \
	--global_param ${DATA}/stratas_input_files/${GLOBAL} \
	--local_param ${DATA}/CNV/${LOCAL} \
	--peaks ${GROUP}/LNCaP/knownGene_withAlias_withgencode.chr${CHR}.txt \
	--cell_specific TRUE \
	--cell_pop ${CELL} \
	--pheno FALSE \
	--gene_express ${DATA}/${EXPRESS} \
	--bbreg TRUE \
	--data_path ${DATA}/stratas_results/${OUT_FOLDER} \
	--betabinom ${BINOM} \
	--sample_name ${SAMPLE}.chr${CHR}.${GENE} \
	--indiv ${INDIV} \
	--gene_name ${GENE} #for running a single gene

termin=$(date +"%s")
difftimelps=$(($termin-$begin))
echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed for Script Execution. #Completed: ${set_num}"
