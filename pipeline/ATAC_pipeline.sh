#!/bin/bash
#SBATCH -n 8                # Number of cores (processes/threads?)
#SBATCH -N 1                # Ensure that all cores are on one machine; number of nodes
#SBATCH --mem=70000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o job_readouts/myoutput_%j.out  #"outFile"$1".txt"# File to which STDOUT will be written, %j inserts jobid
#SBATCH -e job_readouts/myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cakalita@gmail.com

#work from home folder
#stored big files in group
GROUP="../../agusevlab"
AHS="${GROUP}/ashetty"
#WASP="${GROUP}/APPS/WASP"
WASP="${GROUP}/ckal/WASP-master"
#DATA="/DATA/PRCA_CHIP" ##INPUT THIS AS VARIABLE IN QSUB, starting from GROUP directory

#FASTQ1=$1 #starting from DATA directory
#FASTQ2=$2
#VCF=$3 #DNA genotype data, starting from GROUP directory
#SAMPLE=$4 #sample name 

#load modules
module load shared R/3.5.2 gcc/7.3.0 
module load bwa/0.7.17-r1188 samtools/1.9.1 bcftools/1.9 bedtools/2.27.1 python/3.6.7
echo ${SAMPLE} processing has started.

mkdir ${DATA}/step_2_output_files
mkdir ${DATA}/step_2_output_files/sort_temp

echo step 2 has started

# STEP 2
#align reads
#if [ ! -s "${DATA}/step_2_output_files/${SAMPLE}.bam" ]; then
echo starting step 2
bwa mem -t 8 ${GROUP}/LNCaP/ucsc.hg19.fasta \
${DATA}/${FASTQ1} ${DATA}/${FASTQ2} | samtools view -S -b -q 10 - > ${DATA}/step_2_output_files/${SAMPLE}.bam
rm ${DATA}/step_2_output_files/${SAMPLE}.sam
#; fi  #output .sai when bwa aln and single end

#I believe this is not used for paired end reads
#/apps/bwa-0.7.16a/bwa samse \
#$AHS/new_reference_files/ucsc.hg19.fasta \
#${DATA}/step_2_output_files/${SAMPLE}.sai \
#${DATA}/${FASTQ1} | samtools view -b -q 10 - > ${DATA}/step_2_output_files/${SAMPLE}.bam

#sort and index
samtools sort -T ${DATA}/step_2_output_files/sort_temp/${SAMPLE} -O bam ${DATA}/step_2_output_files/${SAMPLE}.bam > ${DATA}/step_2_output_files/${SAMPLE}.sorted.bam
samtools index ${DATA}/step_2_output_files/${SAMPLE}.sorted.bam

#else
#echo step 2 bypassed; fi

echo step 3 has started
# STEP 3
mkdir ${DATA}/step_3_output_files
#make comma delimited list of samples in genotype vcf
bcftools query -l ${GROUP}/${VCF} > ${DATA}/step_3_output_files/lncap_samples.txt
#echo "combined_snp_filtrated" >  ${DATA}/step_3_output_files/combined_snp_filtrated.txt

echo ${DATA}/step_3_output_files/lncap_samples.txt

#wasp run to remove ref bias
if [ ! -s "${DATA}/step_3_output_files/${SAMPLE}.sorted.remap.fq1.gz" ]; then
python3 $WASP/mapping/find_intersecting_snps.py \
	--is_paired_end \
	--is_sorted \
	--output_dir ${DATA}/step_3_output_files \
	--snp_tab ${GROUP}/LNCaP/h5_splitchr_nodup/snp_tab.h5 \
	--snp_index ${GROUP}/LNCaP/h5_splitchr_nodup/snp_index.h5 \
	--haplotype ${GROUP}/LNCaP/h5_splitchr_nodup/haplotypes.h5 \
	--samples ${DATA}/step_3_output_files/lncap_samples.txt \
	${DATA}/step_2_output_files/${SAMPLE}.sorted.bam
echo wasp finished
else
echo wasp failed; fi

echo step 4 has started
# STEP 4
mkdir ${DATA}/step_4_output_files
mkdir ${DATA}/step_4_output_files/sort_temp

#remap using wasp input
if [ ! -s "${DATA}/step_4_output_files/${SAMPLE}.sorted.bam" ]; then
bwa mem -t 8 ${GROUP}/LNCaP/ucsc.hg19.fasta \
${DATA}/step_3_output_files/${SAMPLE}.sorted.remap.fq1.gz \
${DATA}/step_3_output_files/${SAMPLE}.sorted.remap.fq2.gz | samtools view -S -b -q 10 - > ${DATA}/step_4_output_files/${SAMPLE}.bam #.sai for aln

#sort and index
samtools sort -T ${DATA}/step_4_output_files/sort_temp/${SAMPLE} -O bam ${DATA}/step_4_output_files/${SAMPLE}.bam > ${DATA}/step_4_output_files/${SAMPLE}.sorted.bam
samtools index ${DATA}/step_4_output_files/${SAMPLE}.sorted.bam

echo finished realignment
else
echo failed realignment; fi

echo step 5 has started
# STEP 5
mkdir ${DATA}/step_5_output_files
mkdir ${DATA}/step_5_output_files/sort_temp
if [ ! -s "${DATA}/step_5_output_files/${SAMPLE}.keep.bam" ]; then
${GROUP}/APPS/miniconda2/bin/python2.7 $WASP/mapping/filter_remapped_reads.py \
	${DATA}/step_3_output_files/${SAMPLE}.sorted.to.remap.bam \
	${DATA}/step_4_output_files/${SAMPLE}.sorted.bam \
	${DATA}/step_5_output_files/${SAMPLE}.keep.bam
echo finished remapping
else
echo failed remapping; fi

echo step 6 has started
# STEP 6
mkdir ${DATA}/step_6_output_files
#-f indicates existance of file. -s indicates if a file is empty
if [ ! -s "${DATA}/step_6_output_files/${SAMPLE}.keep.merge.bam" ]; then
samtools merge -f ${DATA}/step_6_output_files/${SAMPLE}.keep.merge.bam \
              ${DATA}/step_5_output_files/${SAMPLE}.keep.bam  \
              ${DATA}/step_3_output_files/${SAMPLE}.sorted.keep.bam
echo finished merging for step 6
else
echo failed step 6; fi

echo step 7 has started
# STEP 7
# Be aware the step_6 file is not truly sortted!
mkdir ${DATA}/step_7_output_files
samtools sort -T ${DATA}/step_6_output_files/${SAMPLE} -O bam ${DATA}/step_6_output_files/${SAMPLE}.keep.merge.bam > ${DATA}/step_6_output_files/${SAMPLE}.keep.merge.sorted.bam
samtools index ${DATA}/step_6_output_files/${SAMPLE}.keep.merge.sorted.bam 

if [ -s "${DATA}/step_6_output_files/${SAMPLE}.keep.merge.sorted.bam" ]; then
echo sorting finished
else
echo sorting failed; fi

python3 $WASP/mapping/rmdup_pe.py ${DATA}/step_6_output_files/${SAMPLE}.keep.merge.sorted.bam ${DATA}/step_7_output_files/${SAMPLE}.bam

if [ -s "${DATA}/step_7_output_files/${SAMPLE}.bam" ]; then
echo rmdup finished
else
echo rmdup failed; fi

samtools sort -T ${DATA}/step_7_output_files/${SAMPLE} -O bam ${DATA}/step_7_output_files/${SAMPLE}.bam > ${DATA}/step_7_output_files/${SAMPLE}.sorted.bam
samtools index ${DATA}/step_7_output_files/${SAMPLE}.sorted.bam

if [ -s "${DATA}/step_7_output_files/${SAMPLE}.sorted.bam" ]; then
echo ${DATA}/step_7_output_files/${SAMPLE}.sorted.bam exists and is not empty. finished step1
else
echo failed step7; fi

echo starting ase read counter
mkdir ${DATA}/ase_read_counter_output
java -jar ${GROUP}/ckal/picard.jar AddOrReplaceReadGroups \
I=${DATA}/step_7_output_files/${SAMPLE}.sorted.bam \
O=${DATA}/ase_read_counter_output/${SAMPLE}.red.sorted.bam \
RGLB=lib \
RGPL=illumina \
RGPU=run \
RGSM=${SAMPLE}

rm ${DATA}/ase_read_counter_output/${SAMPLE}.red.sorted.bam.bai
samtools index ${DATA}/ase_read_counter_output/${SAMPLE}.red.sorted.bam

if [ -s "${DATA}/ase_read_counter_output/${SAMPLE}.red.sorted.bam" ]; then
echo ase read counter finished
else
echo ase read counter failed; fi

mkdir ${DATA}/allele_counts
java -Djava.io.tmpdir=/tmp -jar ${GROUP}/ckal/GenomeAnalysisTK-3.8.1/GenomeAnalysisTK.jar \
   -R ${GROUP}/LNCaP/ucsc.hg19.fasta \
   -T ASEReadCounter \
   -o ${DATA}/allele_counts/${SAMPLE}.csv \
   -I ${DATA}/ase_read_counter_output/${SAMPLE}.red.sorted.bam \
   -sites ${GROUP}/${VCF} \
   -U ALLOW_N_CIGAR_READS \
   -minDepth 10

if [ -s "${DATA}/allele_counts/${SAMPLE}.csv" ]; then
echo ${DATA}/allele_counts/${SAMPLE}.csv exists and is not empty. 
else
echo failed to make table of allele counts; fi

echo started making stratas prep files
mkdir ${DATA}/allele_counts
#this grabs the genotypes from the genotype vcf and creates a sorted bed file (header=chr,pos0,pos1,rsid,genotype)
bcftools query -f '%CHROM\t%POS0\t%END\t%ID[\t%GT]\n' ${GROUP}/${VCFname}.vcf \
| sort -k1,1 -k2,2n > ${GROUP}/${VCFname}.bed

#make a sorted bed file from sample with header=chr,pos0,pos1,rsid,refcounts,altcounts
less ${DATA}/allele_counts/${SAMPLE}.csv | awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5,$6,$7}' \
| sed 1d | sort -k1,1 -k2,2n > ${DATA}/allele_counts/${SAMPLE}.bed 

#merge genotypes and sample allele counts
bedtools intersect -a ${GROUP}/${VCFname}.bed -b ${DATA}/allele_counts/${SAMPLE}.bed -wa -wb -sorted > ${DATA}/allele_counts/${SAMPLE}_geno.bed

if [ -s "${DATA}/allele_counts/${SAMPLE}_geno.bed" ]; then
echo intersection of genotypes and allele counts finished
else
echo intersection of genotypes and allele counts failed; fi

mkdir ${DATA}/stratas_prep_files
#remove extra columns, homozygous sites and filter for counts >0. POS is POS1. 1,3-5 are variant data and genotypes, 10-13 are alleles and allele counts. 
less ${DATA}/allele_counts/${SAMPLE}_geno.bed | cut -f1,3-5,10-13 | awk '$4 != "1|1" && $4 != "0|0"' \
| awk -v OFS='\t' 'BEGIN { print "CHR\tPOS\tRSID\tHAP\tREF\tALT\tREF.READS\tALT.READS" } $7 + $8 > 0' \
| awk -v OFS='\t' '{print $1,$2,$3,$5,$6,$4,$7,$8}' > ${DATA}/stratas_prep_files/${SAMPLE}.counts

mkdir ${DATA}/stratas_input_files

Rscript ${GROUP}/APPS/stratAS/params.R \
--min_cov 5 \
--inp_counts ${DATA}/stratas_prep_files/${SAMPLE}.counts \
--out ${DATA}/stratas_input_files/${SAMPLE}

set +o posix

paste <(printf "%s\n%s" ID "${SAMPLE}") <(cat ${DATA}/stratas_input_files/${SAMPLE}.global.params) > ${DATA}/stratas_input_files/${SAMPLE}.modified.global.params

if [ -s "${DATA}/stratas_input_files/${SAMPLE}.global.params" ]; then
echo ${DATA}/stratas_input_files/${SAMPLE}.global.params exists and is not empty. step 3 finished
else
echo stratas param step failed; fi

echo started stratas run

#need matrix of variant info, genotypes, and ref/alt counts, all split by space and then I think it splits lines?
less ${DATA}/stratas_prep_files/${SAMPLE}.counts | sed 's/[|,]/ /g' | sed 1d > ${DATA}/stratas_input_files/${SAMPLE}.OUT.MAT.split
#do the split if you need to parallelize many runs
#| split -d -l 15000 - ${DATA}/stratas_input_files/${SAMPLE}.OUT.MAT.split.

mkdir ${DATA}/stratas_results
#NOTE that peak file needs to be preprocessed to have columns: CHR     P0      P1      NAME    CENTER
Rscript slurm_src/stratas.R \
--input1 ${DATA}/stratas_input_files/${SAMPLE}.OUT.MAT.split \
--global_param ${DATA}/stratas_input_files/${SAMPLE}.modified.global.params \
--peaks ${GROUP}/${PEAK} \
--data_path ${DATA}/stratas_results/ \
--sample_name ${SAMPLE}

#--samples SAMPLES.ID \
#--local_param ${SAMPLE}.local.params \

if [ -s "${DATA}/stratas_results/${SAMPLE}_stratas_results.txt" ]; then
echo ${DATA}/stratas_input_files/${SAMPLE}_stratas_results.txt exists and is not empty. stratas run finished
else 
echo stratas run failed; fi

echo 'completed!'