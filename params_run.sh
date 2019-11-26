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
echo started ${SAMPLE}

mkdir ${DATA}/stratas_prep_files
mkdir ${DATA}/stratas_input_files

Rscript slurm_src/params_new.R \
--min_cov 5 \
--inp_counts ${DATA}/stratas_prep_files/${SAMPLE}.counts \
--out ${DATA}/stratas_input_files/${SAMPLE}
#--inp_cnv ${DATA}/CNV/${LOCAL} \
#--multi_ind TRUE

less ${DATA}/stratas_input_files/${SAMPLE}.global.params | awk -v OFS='\t' '$1 != "NA"' > ${DATA}/stratas_input_files/${SAMPLE}.nona.global.params

termin=$(date +"%s")
difftimelps=$(($termin-$begin))
echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed for Script Execution. #Completed: ${set_num}"
