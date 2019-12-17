#!/bin/bash
#SBATCH -n 8                # Number of cores (processes/threads?)
#SBATCH -N 1                # Ensure that all cores are on one machine; number of nodes
#SBATCH --mem=65000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o job_readouts/myoutput_%j.out  #"outFile"$1".txt"# File to which STDOUT will be written, %j inserts jobid
#SBATCH -e job_readouts/myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cakalita@gmail.com

# This is the STAR align with WASP pipeline script
# Python 2.7, pytables 2.x, pysam and numpy are all required depdendencies for this - install them!
# If my minoconda installation in the lab directory is safe, then these should all be installed already
# This currently requires a folder for the reference fasta - if this does not exist, then this will need to be downloaded
# step{2..7}_output_files/ folders also are needed - hopefully this can be implemented so that we will have only one folder for WASP  

#work from home folder
#stored big files in group
GROUP="../../agusevlab"
WASP="${GROUP}/ckal/WASP-master/"

module load STAR/2.6.1c bwa/0.7.17-r1188 samtools/1.9.1 bcftools/1.9 bedtools/2.27.1 python/3.6.7 R/3.5.1

echo starting ${SAMPLE} processing
#run this once
#apps/STAR-2.6.0a/STAR --runMode genomeGenerate --runThreadN 24 --genomeDir $AHS/new_reference_files/ -genomeFastaFiles $AHS/new_reference_files/ucsc.hg19.fasta
#check this out for doing gene counts https://ucdavis-bioinformatics-training.github.io/2017-June-RNA-Seq-Workshop/wednesday/alignment.html
#using gtf file
#create annotation file
#if [ ! -s "${GROUP}/LNCaP/knownGene.gtf" ]; then
#zcat ${GROUP}/LNCaP/knownGene.txt.gz | cut -f1-10 | ${GROUP}/ckal/genePredToGtf file stdin ${GROUP}/LNCaP/knownGene.gtf; fi

mkdir ${DATA}/STAR
echo ${GROUP}/STAR_index
echo ${DATA}/${FASTQ1}
echo ${DATA}/${FASTQ2}
echo ${GROUP}/${VCF}

STAR --runMode alignReads \
--genomeDir ${GROUP}/STAR_index \
--runThreadN 8 \
--readFilesCommand zcat \
--readFilesIn ${DATA}/${FASTQ1} ${DATA}/${FASTQ2} \
--waspOutputMode SAMtag \
--varVCFfile ${GROUP}/${VCF} \
--outFileNamePrefix ${DATA}/STAR/${SAMPLE}_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI AS nM NM MD vA vG vW

#--sjdbGTFfile ${GROUP}/LNCaP/knownGene.gtf \
if [ -s "${DATA}/STAR/${SAMPLE}_Aligned.sortedByCoord.out.bam" ]; then
echo finished alignment
else
echo failed alignment; fi

samtools view -H ${DATA}/STAR/${SAMPLE}_Aligned.sortedByCoord.out.bam > ${DATA}/STAR/${SAMPLE}.output.sam 
samtools view ${DATA}/STAR/${SAMPLE}_Aligned.sortedByCoord.out.bam | grep 'vW:i:1' >> ${DATA}/STAR/${SAMPLE}.output.sam 
samtools view -S -b -q10 ${DATA}/STAR/${SAMPLE}.output.sam > ${DATA}/STAR/${SAMPLE}_WASP.bam

#sort and index
samtools sort -T ${DATA}/STAR/${SAMPLE} -O bam ${DATA}/STAR/${SAMPLE}_WASP.bam > ${DATA}/STAR/${SAMPLE}_WASP.sorted.bam
samtools index ${DATA}/STAR/${SAMPLE}_WASP.sorted.bam

if [ -s "${DATA}/STAR/${SAMPLE}_WASP.sorted.bam" ]; then
echo step1 finished
else
echo step 1 failed; fi

python3 $WASP/mapping/rmdup_pe.py ${DATA}/STAR/${SAMPLE}_WASP.sorted.bam ${DATA}/STAR/${SAMPLE}_WASP.rmdup.bam

if [ -s "${DATA}/STAR/${SAMPLE}_WASP.rmdup.bam" ]; then
echo rmdup finished
else
echo rmdup failed; fi

#remove sorted bam that didn't have dups removed. kept name since that is used as name in next script
rm ${DATA}/STAR/${SAMPLE}_WASP.sorted.bam
rm ${DATA}/STAR/${SAMPLE}_WASP.sorted.bam.bai
#sort and index
samtools sort -T ${DATA}/STAR/${SAMPLE} -O bam ${DATA}/STAR/${SAMPLE}_WASP.rmdup.bam > ${DATA}/STAR/${SAMPLE}_WASP.sorted.bam
samtools index ${DATA}/STAR/${SAMPLE}_WASP.sorted.bam

mkdir ${DATA}/STAR/ase_read_counter_output
java -jar ${GROUP}/ckal/picard.jar AddOrReplaceReadGroups \
I=${DATA}/STAR/${SAMPLE}_WASP.sorted.bam \
O=${DATA}/STAR/ase_read_counter_output/${SAMPLE}.red.sorted.bam \
RGLB=lib \
RGPL=illumina \
RGPU=run \
RGSM=${SAMPLE}

samtools index ${DATA}/STAR/ase_read_counter_output/${SAMPLE}.red.sorted.bam

if [ -s "${DATA}/STAR/ase_read_counter_output/${SAMPLE}.red.sorted.bam" ]; then
echo picard finished
else
echo picard failed; fi

mkdir ${DATA}/STAR/allele_counts

java -Djava.io.tmpdir=/tmp -jar ${GROUP}/ckal/GenomeAnalysisTK-3.8.1/GenomeAnalysisTK.jar \
   -R ${GROUP}/LNCaP/ucsc.hg19.fasta \
   -T ASEReadCounter \
   -o ${DATA}/STAR/allele_counts/${SAMPLE}.csv \
   -I ${DATA}/STAR/ase_read_counter_output/${SAMPLE}.red.sorted.bam \
   -sites ${GROUP}/${VCF} \
   -U ALLOW_N_CIGAR_READS \
   -minDepth 10

if [ -s "${DATA}/STAR/allele_counts/${SAMPLE}.csv" ]; then
echo step2 finished
else
echo step 2 failed; fi

echo starting step 3

#this grabs the genotypes from the genotype vcf and creates a sorted bed file (header=chr,pos0,pos1,rsid,genotype)
bcftools query -f '%CHROM\t%POS0\t%END\t%ID[\t%GT]\n' ${GROUP}/${VCFname}.vcf \
| sort -k1,1 -k2,2n > ${GROUP}/${VCFname}.bed

#make a sorted bed file from sample with header=chr,pos0,pos1,rsid,refcounts,altcounts
less ${DATA}/STAR/allele_counts/${SAMPLE}.csv | awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5,$6,$7}' \
| sed 1d | sort -k1,1 -k2,2n > ${DATA}/STAR/allele_counts/${SAMPLE}.bed 

#merge genotypes and sample allele counts
bedtools intersect -a ${GROUP}/${VCFname}.bed -b ${DATA}/STAR/allele_counts/${SAMPLE}.bed -wa -wb -sorted > ${DATA}/STAR/allele_counts/${SAMPLE}_geno.bed

if [ -s "${DATA}/STAR/allele_counts/${SAMPLE}_geno.bed" ]; then
echo intersect finished
else
echo intersect failed; fi

mkdir ${DATA}/STAR/stratas_prep_files
#remove extra columns, homozygous sites and filter for counts >0. POS is POS1. 1,3-5 are variant data and genotypes, 10-13 are alleles and allele counts. 
less ${DATA}/STAR/allele_counts/${SAMPLE}_geno.bed | cut -f1,3-5,10-13 | awk '$4 != "1|1" && $4 != "0|0"' \
| awk -v OFS='\t' 'BEGIN { print "CHR\tPOS\tRSID\tHAP\tREF\tALT\tREF.READS\tALT.READS" } $7 + $8 > 0' \
| awk -v OFS='\t' '{print $1,$2,$3,$5,$6,$4,$7,$8}' > ${DATA}/STAR/stratas_prep_files/${SAMPLE}.counts

mkdir ${DATA}/STAR/stratas_input_files

Rscript ${GROUP}/APPS/stratAS/params.R \
--min_cov 5 \
--inp_counts ${DATA}/STAR/stratas_prep_files/${SAMPLE}.counts \
--out ${DATA}/STAR/stratas_input_files/${SAMPLE}

set +o posix

paste <(printf "%s\n%s" ID "${SAMPLE}") <(cat ${DATA}/STAR/stratas_input_files/${SAMPLE}.global.params) > ${DATA}/STAR/stratas_input_files/${SAMPLE}.modified.global.params

if [ -s "${DATA}/STAR/stratas_input_files/${SAMPLE}.global.params" ]; then
echo step 3 finished
else
echo step 3 failed; fi

echo started stratas run

#need matrix of variant info, genotypes, and ref/alt counts, all split by space and then I think it splits lines?
less ${DATA}/STAR/stratas_prep_files/${SAMPLE}.counts | sed 's/[|,]/ /g' | sed 1d > ${DATA}/STAR/stratas_input_files/${SAMPLE}.OUT.MAT.split
#do the split if you need to parallelize many runs
#| split -d -l 15000 - ${DATA}/stratas_input_files/${SAMPLE}.OUT.MAT.split.

mkdir ${DATA}/STAR/stratas_results
#NOTE that peak file needs to be preprocessed to have columns: CHR     P0      P1      NAME    CENTER
Rscript slurm_src/stratas.R \
--input1 ${DATA}/STAR/stratas_input_files/${SAMPLE}.OUT.MAT.split \
--global_param ${DATA}/STAR/stratas_input_files/${SAMPLE}.modified.global.params \
--peaks ${GROUP}/${PEAK} \
--data_path ${DATA}/STAR/stratas_results/ \
--sample_name ${SAMPLE}

#--samples SAMPLES.ID \
#--local_param ${SAMPLE}.local.params \

if [ -s "${DATA}/STAR/stratas_results/${SAMPLE}_stratas_results.txt" ]; then
echo ${DATA}/STAR/stratas_input_files/${SAMPLE}_stratas_results.txt exists and is not empty. stratas run finished
else
echo stratas run failed; fi

echo completed!
