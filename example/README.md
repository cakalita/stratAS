```
#just run for DPF3
Rscript stratas_DeCAF.R \
--input example/KIRC.ALL.AS.chr14.OUT.MAT.split \
--global_param example/KIRC.ALL.AS.global.params \
--local_param example/KIRC.ALL.AS.CNVLOCAL \
--peaks example/knownGene_withAlias.stratas.txt \
--cell_specific TRUE \
--cell_pop example/KIRC_TIMER.txt \
--pheno FALSE \
--gene_express example/KIRC.rnaseq.fixedsamples.txt \
--bbreg TRUE \
--data_path example/stratas_results/ \
--sample_name KIRC.ALL.AS.chr14.DPF3 \
--betabinom TRUE \
--gene_name DPF3

