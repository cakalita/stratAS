```
#just run for DPF3
Rscript stratas_DeCAF.R \
--input example/mat_DPF3.OUT.MAT.SPLIT \
--global_param example/KIRC.ALL.AS.01A.nona.global.params \
--local_param example/KIRC.ALL.AS.CNVLOCAL \
--peaks example/knownGene_withAlias_withgencode.chr14.txt \
--cell_specific TRUE \
--cell_pop example/KIRC_TIMER.txt \
--pheno FALSE \
--gene_express example/gene_express_DPF3.txt \
--bbreg TRUE \
--data_path example/stratas_results/ \
--sample_name KIRC.ALL.AS.chr14.DPF3 \
--betabinom TRUE \
--gene_name DPF3

