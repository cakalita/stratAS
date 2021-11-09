#!/bin/bash
#SBATCH -n 6                # Number of cores (processes/threads?)
#SBATCH --mem=5000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cakalita@gmail.com

begin=$(date +"%s")

#mem used is less if not running all tests
#old mem set to 12000
#old outputting
# -o job_readouts/myoutput_%j.out  #"outFile"$1".txt"# File to which STDOUT will be written, %j inserts jobid
# -e job_readouts/myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid

#work from home folder
#stored big files in group
GROUP="/agusevlab"

#load modules
module load R/4.0 
module load htslib

echo started ${GENE_SET}
echo ${OUT_FOLDER}

mkdir ${DATA}/stratas_results
mkdir ${DATA}/stratas_results/${OUT_FOLDER}

#NOTE that peak file needs to be preprocessed to have columns: CHR     P0      P1      NAME    CENTER

set_num=$(basename ${GENE_SET} | rev | cut -d"." -f1-2 | rev)
CHR=$(basename ${GENE_SET} | rev | cut -d"." -f2 | rev )

echo $set_num
echo $CHR

if [ ${SAMPLE} == "BRCA" ]; then
OUT="${DATA}/stratas_input_files/${SAMPLE}.${set_num}"
else
OUT="${DATA}/stratas_input_files/small/${SAMPLE}.${set_num}"
fi
echo $OUT 

if [ ! -s "$OUT.mat" ]; then
	echo $OUT.mat is empty
	exit
else

echo input is $OUT.mat
echo peaks is ${GENE_SET}
echo data path is ${DATA}/stratas_results/${OUT_FOLDER}

#if [ -z ${RERUN+x} ]; then RERUN="FALSE"; else echo "var is set to '$RERUN'"; fi
if [ "$RERUN" = TRUE ] && [ "${CELLSPECIFIC}" = TRUE ] && [ ! -s "${DATA}/stratas_results/${OUT_FOLDER}${SAMPLE}.${set_num}.cfQTL.bbreg_stratas_results.txt" ]; then
RERUN="FALSE"
fi


NUM_CPU_CORES=$(nproc --all)
memlimit=$(ulimit -m)
memlimit_adj=$(($memlimit-1000000))
ulimit -u $(($NUM_CPU_CORES/4 * 3)) #Limit "X" process to 50% CPU usage
#ulimit -v $memlimit_adj #Limit "X" process to 50% CPU usage
ulimit -m $memlimit_adj
ulimit -a


Rscript slurm_src/decaf_nocnv.R \
	--input $OUT.mat \
	--global_param ${DATA}/stratas_input_files/${GLOBAL} \
	--local_param ${DATA}/CNV/${LOCAL} \
	--peaks ${GENE_SET} \
	--cell_specific ${CELLSPECIFIC} \
	--cell_pop ${CELL} \
	--pheno FALSE \
	--gene_express ${DATA}/${EXPRESS} \
	--bbreg ${BBREG} \
	--eqtl_vanilla ${EQTL_VANILLA} \
	--bbreg_vanilla ${EQTL_VANILLA} \
	--dopurity ${DOPURITY} \
	--doprs ${DOPRS} \
	--prs ${PRS} \
	--dointeraction ${DOINTERACTION} \
	--purity ${DATA}/${PURITY} \
	--data_path ${DATA}/stratas_results/${OUT_FOLDER} \
	--betabinom ${BINOM} \
	--indiv ${INDIV} \
	--rerun ${RERUN} \
	--ranknorm FALSE \
	--window ${WINDOW} \
	--sample_name ${SAMPLE}.${set_num}
#	--gene_name ${GENE} \

#	--ranknorm TRUE #default is set to false 


termin=$(date +"%s")
difftimelps=$(($termin-$begin))
echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed for Script Execution. #Completed: ${set_num}"

fi
