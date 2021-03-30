"""
Revisit 5 datasets
  1. Prepare mean_psi (nrows = sample) and expr_file (sample * gene name)
  2. Calculate correlation between gene expression and 2 mean psi
  3. GO analysis
  4. the effect of downsampling on 2 mean psi
  5. ssGSEA: find pathway activity correlated to 2 psis (increase the power of the test)
"""
########## 2. human HSC #########
# 1.1 Prepare mean_psi (nrows = sample) and expr_file (sample * gene name)
mean_psi_file = "HSC/HSC_psi_bimodal.csv"
mean_psi <- read.csv(mean_psi_file,header = TRUE,stringsAsFactors = FALSE, sep=" ") 
expr_file <- read.csv("HSC/HSC_normed_filtered.csv",header = TRUE,stringsAsFactors = FALSE, sep=" ") 
expr_file<-as.data.frame(t(expr_file))

# 1.2 Calculate correlation between gene expression and 2 mean psi 
# colnames 中35965个基因有2453个没有完成enID到gene names的转换
# this is one cluster of cells
library(ggplot2)
expr_psi <- merge_by_rownames(expr_file, mean_psi)
e <- find_gene_corr_FDR(expr_psi,"glbl_mean_psi",cor_method="pearson", tail_col=2)
posi_corr_list_data2 <- e[[1]]
nega_corr_list_data2 <- e[[2]]
# several hundreds

e <- find_gene_corr_FDR(expr_psi,"gene_mean_psi",cor_method="pearson", tail_col=2)
posi_corr_list_data2_psi2 <- e[[1]]
nega_corr_list_data2_psi2 <- e[[2]]
# 20-100

# filter out genes that are not DEG
plot_venn(DEG ,nega_corr_list_data2, nega_corr_list_data2_psi2, "HSC/HSC_nega_corr.tiff")

# 1.3 GO analysis
plots1<-GO_analysis(setdiff(posi_corr_list_data2, DEG), colnames(expr_psi))
plots2<-GO_analysis(setdiff(nega_corr_list_data2, DEG), colnames(expr_psi))

# 1.5 ssGSEA: find pathway activity correlated to 2 psis (increase the power of the test)
#BiocManager::install("GSVA")
library(GSVA)
library(GSEABase)
expr_mat=t(expr_file)  # a matrix of expression values where rows correspond to genes and columns correspond to samples.
# .gmt file通过write_geneset_to_csv.ipynb转成.csv，一列一个GO，列名为GO名
# 必须加na.strings = ""， 空值才会被识别为NA
gs = read.csv("C5_go.csv", stringsAsFactors = FALSE, check.names = FALSE,na.strings = "")
gs2 = as.list(gs)
gs3 = lapply(gs2, function(x) x[!is.na(x)])

# limit the size of geneset to increase power
ssgsea_score = gsva(expr_mat, gs3[2:length(gs3)], method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
#ssgsea_score_limit = gsva(expr_mat, gs_limit, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'

write.csv(ssgsea_score, "HSC/HSC_ssGSEA.csv")

ssgsea_score<-t(ssgsea_score)
gsea_psi_bi <- merge_by_rownames(ssgsea_score, mean_psi)

e <- find_gene_corr_FDR(gsea_psi_bi,"glbl_mean_psi",cor_method="pearson",cor_threshold=0.3,tail_col=2)
posi_gsea_data2_bi <- e[[1]]
nega_gsea_data2_bi <- e[[2]]

e <- find_gene_corr_FDR(gsea_psi_bi,"gene_mean_psi",cor_method="pearson",tail_col=2)
posi_gsea_data2_bi_psi2<- e[[1]]
nega_gsea_data2_bi_psi2 <- e[[2]]
# None of the pathway is correlated to both

plot_venn(DEG ,posi_gsea_data2_bi ,posi_gsea_data2_bi_psi2,"HSC/HSC_posi_corr_gsea.tiff")

e<- find_gene_corr_FDR(gsea_psi_bi,"glbl_mean_psi",cor_method="pearson",tail_col=2)




