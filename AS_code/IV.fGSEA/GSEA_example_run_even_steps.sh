#!/bin/sh

cat GSEA_fifth.txt | parallel eval {}


#python GSEA_class.py mouseSC2/mouseSC2_qNSC1_GENIE3_rank.csv c5.all.v7.2.symbols.gmt mouseSC2/mouseSC2_qNSC1_GENIE3_gene_set_list 100000

#python GSEA_class.py mouseSC2/mouseSC2_qNSC1_corr_rank.csv c5.all.v7.2.symbols.gmt mouseSC2/mouseSC2_qNSC1_corr_positive_gene_set_list 100000

#python GSEA_class.py mouseSC2/mouseSC2_qNSC1_corr_rank.csv c5.all.v7.2.symbols.gmt mouseSC2/mouseSC2_qNSC1_corr_negative_gene_set_list 100000 -inv 1


#python GSEA_class.py mouseSC/mouseSC_18W_GENIE3_rank.csv c5.all.v7.2.symbols.gmt mouseSC/mouseSC_18W_GENIE3_gene_set_list 100000

#python GSEA_class.py mouseSC/mouseSC_18W_corr_rank.csv c5.all.v7.2.symbols.gmt mouseSC/mouseSC_18W_corr_positive_gene_set_list 100000

#python GSEA_class.py mouseSC/mouseSC_18W_corr_rank.csv c5.all.v7.2.symbols.gmt mouseSC/mouseSC_18W_corr_negative_gene_set_list 100000 -inv 1

#python GSEA_class.py iPSC/NPC_GENIE3_rank_bimodal.csv c5.all.v7.2.symbols.gmt iPSC/NPC_GENIE3_gene_set_list 100000

#python GSEA_class.py iPSC/NPC_corr_rank_bimodal.csv c5.all.v7.2.symbols.gmt iPSC/NPC_corr_positive_gene_set_list 100000

#python GSEA_class.py iPSC/NPC_corr_rank_bimodal.csv c5.all.v7.2.symbols.gmt iPSC/NPC_corr_negative_gene_set_list 100000 -inv 1

#python GSEA_class.py iPSC/MN_GENIE3_rank_bimodal.csv c5.all.v7.2.symbols.gmt iPSC/MN_GENIE3_gene_set_list 100000

#python GSEA_class.py iPSC/MN_corr_rank_bimodal.csv c5.all.v7.2.symbols.gmt iPSC/MN_corr_positive_gene_set_list 100000

#python GSEA_class.py iPSC/MN_corr_rank_bimodal.csv c5.all.v7.2.symbols.gmt iPSC/MN_corr_negative_gene_set_list 100000 -inv 1


#python GSEA_class.py hESC/hESC_primed_GENIE3_rank.csv c5.all.v7.2.symbols.gmt hESC/hESCp_GENIE3_gene_set_list 100000

#python GSEA_class.py hESC/hESCp_corr_rank.csv c5.all.v7.2.symbols.gmt hESC/hESCp_corr_positive_gene_set_list 100000

#python GSEA_class.py hESC/hESCp_corr_rank.csv c5.all.v7.2.symbols.gmt hESC/hESCp_corr_negative_gene_set_list 100000 -inv 1

#python GSEA_class.py hESC/hESC_naive_GENIE3_rank.csv c5.all.v7.2.symbols.gmt hESC/hESCn_GENIE3_gene_set_list 100000

#python GSEA_class.py hESC/hESCn_corr_rank.csv c5.all.v7.2.symbols.gmt hESC/hESCn_corr_positive_gene_set_list 100000

#python GSEA_class.py hESC/hESCn_corr_rank.csv c5.all.v7.2.symbols.gmt hESC/hESCn_corr_negative_gene_set_list 100000 -inv 1

#python GSEA_class.py mouseSC/mouseSC_GENIE3_rank.csv c5.all.v7.2.symbols.gmt mouseSC/mouseSC_GENIE3_gene_set_list 100000

#python GSEA_class.py mouseSC/mouseSC_corr_rank.csv c5.all.v7.2.symbols.gmt mouseSC/mouseSC_corr_positive_gene_set_list 100000

#python GSEA_class.py mouseSC/mouseSC_corr_rank.csv c5.all.v7.2.symbols.gmt mouseSC/mouseSC_corr_negative_gene_set_list 100000 -inv 1

#python GSEA_class.py mouseSC2/mouseSC2_GENIE3_rank.csv c5.all.v7.2.symbols.gmt mouseSC2/mouseSC2_GENIE3_gene_set_list 100000

#python GSEA_class.py mouseSC2/mouseSC2_corr_rank.csv c5.all.v7.2.symbols.gmt mouseSC2/mouseSC2_corr_positive_gene_set_list 100000

#python GSEA_class.py mouseSC2/mouseSC2_corr_rank.csv c5.all.v7.2.symbols.gmt mouseSC2/mouseSC2_corr_negative_gene_set_list 100000 -inv 1

#python GSEA_class.py iPSC/iPSC_GENIE3_rank_bimodal.csv c5.all.v7.2.symbols.gmt iPSC/iPSC_GENIE3_gene_set_list 100000

#python GSEA_class.py iPSC/iPSC_corr_rank_bimodal.csv c5.all.v7.2.symbols.gmt iPSC/iPSC_corr_positive_gene_set_list 100000

#python GSEA_class.py iPSC/iPSC_corr_rank_bimodal.csv c5.all.v7.2.symbols.gmt iPSC/iPSC_corr_negative_gene_set_list 100000 -inv 1

#python GSEA_class.py HSC/HSC_GENIE3_rank_bimodal.csv c5.all.v7.2.symbols.gmt HSC/HSC_GENIE3_gene_set_list 100000

#python GSEA_class.py HSC/HSC_corr_rank_bimodal.csv  c5.all.v7.2.symbols.gmt HSC/HSC_corr_positive_gene_set_list 100000 

#python GSEA_class.py HSC/HSC_corr_rank_bimodal.csv  c5.all.v7.2.symbols.gmt HSC/HSC_corr_negative_gene_set_list 100000 -inv 1~                                                                                                                                      


# Give 2 arguments, input and utput fi]es for FDR correction  
#Rscript FDR_correct.R iPSC_corr_negative_gene_set_list_20000.p.csv iPSC_corr_negative_20000_FDR.p.csv
#Rscript FDR_correct.R HSC_corr_negative_20000_gene_set_list.txt.p_value.csv HSC_corr_negative_20000_FDR_p_value.csv 
