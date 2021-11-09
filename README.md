# :first_quarter_moon: DeCAF(DEconvoluted Cell type Allele specific Function)

Identify cell-type specific QTL effects in tumors by leveraging both allelic and total expression information

## Workflow

### Input data / dependencies

* A `.vcf` file containing phased individual genotypes and an "AS" field listing the allelic counts. This gets turned into a flat file to be loaded in the `--input` flag (see below).
* When --cell_specific flag is used: A cell fraction file with row per individual (needs to be labeled as sample_ID, with sample_type as additional column), columns per cell type. Also need a gene_express file containing gene expression with "gene name" and then columns for expression in each individual (labeled).
* When --purity flag is used: A tumor purity file with row per individual (needs to be labeled as sample_ID, with sample_type as additional column), columns per cell type. ie/CPE estimates from Aran et al 2015 10.1038/ncomms9971
* When --pheno flag is used: A samples file with the columns `ID` and `CONDITION` that is in the **same order** as the vcf input. `CONDITION` is currently restricted to 0/1. This file got into the `--samples` flag. This is currently not supported in combination with other conditions.
* A `.bed` file listing the physical positions of the features to be tested for allelic imbalance. This goes into the `--peaks` flag. ie/ gene positions
* Global parameter files generated by `params_new.R` (see below). At minimum a global overdispersion parameter file with the columns "ID","CHR","P0","P1","PHI","MU","N" with entries in any order goes into the `--global_params` flag.
* Local parameter (CNV) generated by `PARAMS_CNV.R` goes into `local_param` flag.
* R, and the `VGAM`, `optparse`, `plyr`, `parallel`, `dplyr`, and `data.table` libraries installed.

### Recommended pre-processing and allelic counts

* Phase your genotype data with EAGLE/HRC.
* Process all sequence data with the [WASP](https://github.com/bmvdgeijn/WASP) mapping pipeline.
* Use the GATK ASEReadCounter to count reads then convert into vcf (see `pipeline/` scripts).

### `params_new.R` : estimating prior parameters

`params_new.R` infers the read count distribution (beta binomial) parameters and needs to be run on a whole genome vcf one individual at a time. The input is an `--inp_counts` file which must contain the columns - 1:4= chr,pos,RSID,ref_allele then every 4 are sample,hap,ref.reads,alt.reads where `HAP` is the `0|1` or `1|0` vcf haplotype code.

A vcf file containing one individual is converted to counts as follows:
```
bcftools query -f '%CHROM\t%END\t%ID\t%REF[\t%SAMPLE\t%GT\t%AS]\n' KIRC.ALL.AS.filtered_01A.vcf.gz | tr ',' '\t' > KIRC.ALL.AS.01A.counts
```

For tumor data, local CNV-specific parameters are also estimated, and the `--inp_cnv` file must be provided with headers `CHR P0 P1` listing the boundaries of CNV regions. This file can additionally contain a `CNV` column listing the CNV estimate (centered to zero as in TCGA calls) for inclusion as a covariate in the final analysis.

inference is then performed by running:
```
sbatch --job-name=params --export=SAMPLE="KIRC.ALL.AS.01A",DATA="../../agusevlab/ckal/TCGA_vcf/KIRC" slurm_src/stratas_params_CGC.sh 
```
Where the content looks like:
```
Rscript slurm_src/params_new.R \
--min_cov 5 \
--inp_counts ${DATA}/stratas_prep_files/${SAMPLE}.counts \
--out ${DATA}/stratas_input_files/${SAMPLE}
#--inp_cnv ${DATA}/CNV/${LOCAL} \
```
Uncomment --inp_cnv when calculating local CNV (this is replaced by the `PARAMS_CNV.R` script). Then remove NA from the params file which can cause issues in the stratas run.
```
less ${DATA}/stratas_input_files/${SAMPLE}.global.params | awk -v OFS='\t' '$1 != "NA"' > ${DATA}/stratas_input_files/${SAMPLE}.nona.global.params

```

A `$OUT.global.params` file is generated containing the parameters (with header `PHI MU N`) and (optionally) an `$OUT.local.params` file is generated containing the positions and parameters for each CNV (with header `ID CHR P0 P1 PHI MU N`).

### `decaf_nocnv.R` : testing for differential AS

`decaf_nocnv.R` computes the AS statistics from a VCF of all individuals and the prior parameters estimated above.

To split gene/peak file into chromosomes, and then into smaller sets:
```
awk -v OFS='\t' '{print > ("knownGene_withAlias_withgencode."$1".txt")}' knownGene_withAlias_withgencode.txt

for chr in {1..22}; do 
split -d -l 10 -a 3 knownGene_withAlias_withgencode.chr${chr}.txt knownGene_withAlias_withgencode.chr${chr}.
done
```

A vcf file containing all individuals is converted to counts and split into chromosomes as follows:
```
bcftools query -f '%CHROM\t%END\t%ID\t%REF\t%ALT[\t%GT\t%AS]\n' KIRC.ALL.AS.filtered_01A.vcf.gz | tr ',' '\t' | sed 's/[|,]/\t/g' | sed 1d > ../stratas_input_files/KIRC.ALL.AS.01A.OUT.MAT.split
awk -v OFS='\t' -v var="KIRC.ALL.AS.01A" '{print > (var"."$1".OUT.MAT.split")}' KIRC.ALL.AS.01A.OUT.MAT.split
```

Each batch is then processed as follows:
```
for gene_set in ../../agusevlab/LNCaP/split_gene/*chr*+([0-9]); do
  sbatch --job-name=${gene_set} --export=BINOM="TRUE",LOCAL="KIRC.ALL.AS.CNVLOCAL",GLOBAL="KIRC.ALL.AS.nona.global.params",SAMPLE="KIRC.ALL.AS.01A",DATA="../../agusevlab/ckal/TCGA_vcf/KIRC",PEAK="LNCaP/split_gene/${gene_set}",GENE_SET="${gene_set}",CELL="../../agusevlab/ckal/AIsim/KIRC_TIMER.txt",OUT="batch2/",EXPRESS="broad_rnaseq/KIRC.rnaseq.fixedsamples.full.txt" slurm_src/batch_decaf.sh 
  sleep 1
done
```
### Output data

By default, the ASE test is printed to file testtype_stratas_results.txt with each line containing the following entries:

| Column | Description |
| --- | --- |
| CHR | Chromosome |
| POS | Position of test SNP |
| RSID | ID of test SNP |
| P0 | Start of gene/peak  |
| P1 | End of gene/peak |
| NAME | Name of gene/peak |
| CENTER | Center position of peak (or TSS for gene) |
| N.HET | # of heterozygous individuals tested |
| N.READS | # of reads tested in total |

Then additional columns depending on the test type choice/s made. 
For -betabinom
| Column | Description |
| --- | --- |
| ALL.AF | Allelic fraction estimate from standard beta binomial test |
| ALL.BBINOM.P | Beta-binomial test for imbalance across both conditions  |
| SAMPLE_NAME | Name of the Sample |

For -bbreg 
| Column | Description |
| --- | --- |
| z | z score from beta binomial test that includes CNV as a covariate |
| pv | p value from beta binomial test that includes CNV as a covariate |
| SAMPLE_NAME | Name of the Sample |

For -cell_specific you get 2 files. eqtl results were added for DeCAF, so they are currently only output when this flag is true although cell pop data is not input.

eqtl_vanilla results:
| Column | Description |
| --- | --- |
| z_eqtl_vanilla | z score from eqtl interaction test that includes CNV |
| eqtl_pval_vanilla | p value from eqtl interaction test that includes CNV |
| SAMPLE_NAME | Name of the Sample |

DeCAF results (this output is the same for purity and cell fraction calculations):
| Column | Description |
| --- | --- |
| z_eqtl | z score from eqtl interaction test that includes CNV and cell fraction |
| z_eqtl_pval | p value from eqtl interaction test that includes CNV and cell fraction |
| z_AI | z score from bbreg betabinomial interaction test that includes CNV and cell fraction |
| z_AI_pval | p value from bbreg betabinomial interaction test that includes CNV and cell fraction |
| c | number of tests with non-NA results (used for Stouffer's combination calculation) |
| z_comb | z score from combining eqtl and bbreg betabinomial tests with Stouffer's method |
| z_comb_pval | p value from combining eqtl and bbreg betabinomial tests with Stouffer's method |
| sample | Name of the Sample |

For -pheno and -betabinom
| Column | Description |
| --- | --- |
| ALL.AF | Allelic fraction estimate from beta binomial test across both conditions |
| ALL.BBINOM.P | Beta-binomial test for imbalance across both conditions  |
| C0.AF | Allelic fraction estimate from condition 0 |
| C0.BBINOM.P | Beta-binomial test for imbalance in condition 0 |
| C1.AF | Allelic fraction estimate from condition 1 |
| C1.BBINOM.P | Beta-binomial test for imbalance in condition 1 |
| DIFF.BBINOM.P | Beta-binomial test for difference between conditions |

Enabling the `--binom` flag additionally runs a standard binomial test, and produces the following columns:
| Column | Description |
| --- | --- |
| ALL.BINOM.P | Binomial test for imbalance across both conditions|
| ALL.C0.BINOM.P | Binomial test for imbalance in condition 0 |
| ALL.C1.BINOM.P | Binomial test for imbalance in condition 1 |
| FISHER.OR | Fisher's test odd's ratio for difference between conditions|
| FISHER.DIFF.P | Fisher's test difference between conditions |

Enabling the -pheno and --bbreg flag runs a beta binomial regression with CNV as covariate, and produces the following columns:
| Column | Description |
| --- | --- |
| ALL.BBREG.P | Beta binomial regression (with covariates) for imbalance across both conditions |
| DIFF.BBREG.P | Beta binomial regression (with covariates) for imbalance difference between conditions |
| CNV.BBREG.P | Beta binomial regression (with covariates) for imbalance along CNV covariate |

Enabling the `--indiv` flag additionally produces the following columns:
| Column | Description |
| --- | --- |
| IND.C0 | Number of each condition 0 individual included in this test (comma separated) |
| IND.C0.COUNT.REF | condition 0 REF allele counts of each individual included in this test (comma separated) |
| IND.C0.COUNT.ALT | condition 0 ALT allele counts of each individual included in this test (comma separated) |
| IND.C1 | Number of each condition 1 individual included in this test (comma separated) |
| IND.C1.COUNT.REF | condition 1 REF allele counts of each individual included in this test (comma separated) |
| IND.C1.COUNT.ALT | condition 1 ALT allele counts of each individual included in this test (comma separated) |

## Example

An example locus with significant AS associations can be run by calling:

```
#just run for DPF3
Rscript decaf_nocnv.R \
--input example/KIRC.ALL.AS.01A.chr14.194.mat \
--global_param example/KIRC.ALL.AS.01A.nona.new.global.params \
--local_param example/KIRC.RNA.new.local.params \
--peaks example/knownGene_withAlias_withgencode.chr14.194.DPF3 \
--cell_specific TRUE \
--cell_pop example/KIRC_TIMER_new.txt \
--pheno FALSE \
--gene_express example/KIRC.rnaseq.fixedsamples.fullID.txt \
--bbreg TRUE \
--eqtl_vanilla TRUE \
--bbreg_vanilla TRUE \
--dopurity TRUE \
--dointeraction TRUE \
--purity example/tumor_purity.txt \
--data_path example/stratas_results/ \
--sample_name KIRC.ALL.AS.chr14.DPF3 \
--betabinom FALSE 

```
## Detailed Parameters

| Parameter | Description |
| --- | --- |
| `--input` | Path to input file |
| `--samples` | Path to sample identifier file, must have ID and CONDITION columns |
| `--peaks` | Path to file containing peak/gene boundaries, must contain CHR P0 P1 NAME CENTER columns |
| `--global_param` | Path to global overdispersion parameter file |
| `--local_param` | Path to local overdispersion parameter file |
| `--out` | Path to output |
| `--window` | Window (in bp) for SNPs to test around the peak boundary |
| `--perm` | # of permutations to shuffle the allele labels (0=off) |
| `--perm_cond` | # of permutations to shuffle the condition labels (0=off) |
| `--min_cov` | Individuals must have at least this many reads (for both alleles) to be tested |
| `--min_maf` | Minimum minor allele frequency for test SNP |
| `--min_het` | Minimum minor heterozygous frequency for test SNP |
| `--max_rho` | Maximum local/global over-dispersion parameter for which to include individual in test |
| `--binom` | Also perform a standard binomial test |
| `--bbreg` | Also perform a beta binomial regression with local CNV status as a covariate (must also provide `--local_param` file) |
| `--indiv` | Also report the per-individual allele fractions (Warning, this can produce large files) |
| `--exclude` | The mimium distance between SNPs allowed in the haplotype (to exclude variants in the same read) |
| `--gene_express` | Path to input file containing gene expression results needed for cell-type specific QTL testing |
| `--sample_name` | Name of library being analyzed |
| `--gene_name` | Name of gene being analyzed |
| `--cell_pop` | Path to file containing cell population fractions, row per individual, columns per cell type |
| `--out` | Path to output |
| `--cell_specific` | Perform cell-type specific ASE/QTL test |
| `--data_path` | Set output data path to store results |
| `--pheno` | Preforms liklihood ratio test for chisq results between conditions. Curently binary only |
| `--betabinom` | Perform a standard beta-binomial test without LRT |
| `--dopurity` | Perform tumor specific ASE/QTL test |
| `--purity` | Path to file containing tumor purity estimates, row per individual |
| `--eqtl_vanilla` | Perform marginal QTL test (includes tumor purity covariate)|
| `--bbreg_vanilla` | Perform marginal AI test |
| `--dointeraction` | Perform any of the interaction tests |
| `--rerun` | If script does not run through the full input, this can be set to TRUE and loop will jump to the last run SNP |


## Notes:

* Permutation: The permutation test permutes read counts with respect to alleles. This is a null of counts randomly sampled from the observed count distribution. If a small number of individuals have unusually high read counts and imbalance they will dominate the observed imbalance and appear highly significant even in permutation. In these cases, failing the permutation test may rule out true biological signal due to individual outliers.

* Unlike an eQTL, we can test for AS signal within each individual separately but this is currently not being evaluated. It may be interesting to know what fraction of individual heterozygotes exhibit significant imbalance (for example, through the Storey & Tibshirani pi statistic). Alternatively, individual standard errors can be approximated from each individual test and heterogeneity assessed using standard meta-analysis statistics.

* For a test SNP that is not in the target feature, it may be of interest to report the allelic imbalance for homozygous individuals (that have heterozygous read-carrying variants). If the test SNP explains all of the imbalance at the locus we would expect these homozygous individuals to follow the null. This could also be a way to identify secondary associations.
