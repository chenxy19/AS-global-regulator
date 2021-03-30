"""
Revisit 5 datasets
  1. Prepare mean_psi (nrows = sample) and expr_file (sample * gene name)
  2. Calculate correlation between gene expression and 2 mean psi
  3. GO analysis
  4. the effect of downsampling on 2 mean psi
  5. ssGSEA: find pathway activity correlated to 2 psis (increase the power of the test)
"""
########## 1. human MN, NPC, iPSC #########
# 1.1 Prepare mean_psi (nrows = sample) and expr_file (sample * gene name)
mean_psi_file = "/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/all_mean_psi_corrected0831.csv"
mean_psi <- read.csv(mean_psi_file,header = TRUE,stringsAsFactors = FALSE, sep=" ") 
expr_file = SR_normed@assays$RNA@data
expr_file <- t(as.matrix(expr_file))
rownames(expr_file) <- gsub("TPM_","",rownames(expr_file))

# load mean psi for  bimodal genes
iPSC2p_2meth <- read.csv("/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/iPSC_2p_2meth.csv",header = TRUE,stringsAsFactors = FALSE, sep=" ") 
MN2p_2meth <- read.csv("/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/MN_psi_bimodal.csv",header = TRUE,stringsAsFactors = FALSE, sep=" ") 
NPC2p_2meth <- read.csv("/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/NPC_psi_bimodal.csv",header = TRUE,stringsAsFactors = FALSE, sep=" ") 


# 1.2 Calculate correlation between gene expression and 2 mean psi 
# colnames 中13677个基因有580个左右没有完成enID到gene names的转换
library(ggplot2)
expr_psi <- merge_by_rownames(expr_file, mean_psi)
e <- find_gene_corr_FDR(expr_psi,"glbl_mean_psi",cor_method="pearson")
posi_corr_list_data1_all <- e[[1]]
nega_corr_list_data1_all <- e[[2]]
# this is a set, not ranked； several hundreds

e <- find_gene_corr_FDR(expr_psi,"mean_psi_genes",cor_method="pearson")
posi_corr_list_data1_all_psi2 <- e[[1]]
nega_corr_list_data1_all_psi2 <- e[[2]]
# this is a set, not ranked；<5个

part_expr_psi <- expr_psi[which(expr_psi$cell_type=="iPSC"),]
e <- find_gene_corr_FDR(part_expr_psi,"glbl_mean_psi",cor_method="pearson")
posi_corr_list_data1_iPSC <- e[[1]]
nega_corr_list_data1_iPSC <- e[[2]]
# Bonferroni和BH correction都没有找到显著的

iPSC_marker <- FindMarkers(SR_normed,ident.1 = "iPSC",slot = "data")
NPC_marker <- FindMarkers(SR_normed,ident.1 = "NPC",slot = "data")
MN_marker <- FindMarkers(SR_normed,ident.1 = "MN",slot = "data")
DEG <- union(rownames(MN_marker), rownames(NPC_marker))
DEG <-union(DEG, rownames(iPSC_marker))

# filter out genes that are not DEG
plot_venn(DEG ,nega_corr_list_data1_all,nega_corr_list_data1_all_psi2)

# 1.2.2 Calculate correlation with bimodal mean psi
MN_expr_psi_bi<-merge_by_rownames(expr_file, MN2p_2meth)
NPC_expr_psi_bi<-merge_by_rownames(expr_file, NPC2p_2meth)
#MN_expr_psi_bi<-MN_expr_psi_bi[,-((ncol(MN_expr_psi_bi)-2):ncol(MN_expr_psi_bi))]

e <- find_gene_corr_FDR(iPSC_expr_psi_bi,"bi_glbl_mean_psi",cor_method="pearson",tail_col=4)
posi_corr_list_data1_bi <- e[[1]]
nega_corr_list_data1_bi <- e[[2]]
# none

e <- find_gene_corr_FDR(iPSC_expr_psi_bi,"bi_gene_mean_psi",cor_method="pearson",tail_col=4)
posi_corr_list_data1_bi_psi2 <- e[[1]]
nega_corr_list_data1_bi_psi2 <- e[[2]]
# none

# write ranking to file
corr_volcano_NPC<-find_gene_high_correlation(NPC_expr_psi_bi,"bi_glbl_mean_psi",
                                     cor_method="pearson",
                                     cor_threshold=0,
                                     tail_col=2)
corr_rank_NPC<-rownames(corr_volcano_NPC)[order(corr_volcano_NPC$corr,decreasing=TRUE)]  #进行排序
write.table (corr_rank_NPC, file = "/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/NPC_corr_rank_bimodal.csv", row.names =TRUE, col.names =TRUE, quote =TRUE)

# 1.3 GO analysis
plots<-GO_analysis(setdiff(posi_corr_list_data1_all, DEG), colnames(expr_psi))
plots<-GO_analysis(setdiff(nega_corr_list_data1_all, DEG), colnames(expr_psi))

# write ranking to file
iPSC_GENIE3_rank <-run_GENIE3(t(iPSC_expr_psi_bi),"bi_glbl_mean_psi","/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/iPSC_GENIE3_rank_bimodal.csv")

MN_expr_psi_bi2<-MN_expr_psi_bi[,-ncol(MN_expr_psi_bi)]
MN_GENIE3_rank <-run_GENIE3(t(MN_expr_psi_bi2),"bi_glbl_mean_psi","/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/MN_GENIE3_rank_bimodal.csv")

NPC_expr_psi_bi2<-NPC_expr_psi_bi[,-ncol(NPC_expr_psi_bi)]
NPC_GENIE3_rank <-run_GENIE3(t(NPC_expr_psi_bi2),"bi_glbl_mean_psi","/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/NPC_GENIE3_rank_bimodal.csv")


# 1.4 the effect of downsampling on 2 mean psi
junc_mat_file="/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/iPSC_junc_mat.csv"
iPSC_junc_mat= read.csv(junc_mat_file,sep=" ",header = TRUE,stringsAsFactors = FALSE) 
iPSC_psi_junc_cell <- calculate_psi("data1_MN+NPC+iPSC/iPSC",iPSC_junc_mat)
library(DropletUtils)
iPSC_junc_mat_down<-downsampleMatrix(as.matrix(iPSC_junc_mat), prop=0.1) #sum of each row *0.5
iPSC_psi_junc_cell_down <- calculate_psi("data1_MN+NPC+iPSC/iPSC",iPSC_junc_mat_down)

# gene-averaged psi
iPSC_gene_mean_psi <- apply(iPSC_psi_junc_cell, 2, mean_non_na)  
iPSC_gene_mean_psi_down <- apply(iPSC_psi_junc_cell_down, 2, mean_non_na)  
cor.test(iPSC_gene_mean_psi, iPSC_gene_mean_psi_down)
# 50%: 0.9396896, 20%: 0.608, 10%: 0.417

# junction_sum psi
iPSC_glbl_mean_psi<-calculate_global_mean_psi("data1_MN+NPC+iPSC/iPSC",iPSC_junc_mat)
iPSC_glbl_mean_psi_down<-calculate_global_mean_psi("data1_MN+NPC+iPSC/iPSC",iPSC_junc_mat_down)
cor.test(iPSC_glbl_mean_psi, iPSC_glbl_mean_psi_down)
# 50%: 0.9960144 , 20%: 0.9226952 , 10%: 0.7557014 

# 1.5 ssGSEA: find pathway activity correlated to 2 psis (increase the power of the test)
#BiocManager::install("GSVA")
library(GSVA)
library(GSEABase)
expr_mat=t(expr_file)  # a matrix of expression values where rows correspond to genes and columns correspond to samples.
# .gmt file通过write_geneset_to_csv.ipynb转成.csv，一列一个GO，列名为GO名
# 必须加na.strings = ""， 空值才会被识别为NA
gs = read.csv("/Users/chenxinyi/Desktop/output/C5_go.csv", stringsAsFactors = FALSE, check.names = FALSE,na.strings = "")
gs2 = as.list(gs)
gs3 = lapply(gs2, function(x) x[!is.na(x)])

# limit the size of geneset to increase power
gs_limit <- gs3[which(as.numeric(lapply(gs3,length)) > 30 & as.numeric(lapply(gs3,length))<200)]
data1_ssgsea_score = gsva(expr_mat, gs3[2:length(gs3)], method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
data1_ssgsea_score_limit = gsva(expr_mat, gs_limit, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'

write.csv(data1_ssgsea_score, "/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/data1_ssGSEA.csv")

data1_ssgsea_score_limit<-t(data1_ssgsea_score_limit)
#gsea_psi <- merge_by_rownames(data1_ssgsea_score_limit, mean_psi)
gsea_psi_bi <- merge_by_rownames(data1_ssgsea_score_limit, iPSC2p_2meth)

e <- find_gene_corr_FDR(gsea_psi_bi,"bi_glbl_mean_psi",cor_method="spearman",tail_col=4)
posi_gsea_data1_bi <- e[[1]]
nega_gsea_data1_bi <- e[[2]]

e <- find_gene_corr_FDR(gsea_psi_bi,"bi_gene_mean_psi",cor_method="spearman",tail_col=4)
posi_gsea_data1_bi_psi2<- e[[1]]
nega_gsea_data1_bi_psi2 <- e[[2]]
# None of the pathway is correlated to both

plot_venn(DEG ,posi_gsea_data1_bi ,posi_gsea_data1_bi_psi2)

part_gsea_psi <- gsea_psi[which(gsea_psi$cell_type=="iPSC"),]
e <- find_gene_corr_FDR(part_gsea_psi,"glbl_mean_psi",cor_method="pearson")
posi_gsea_data1_iPSC <- e[[1]]
nega_gsea_data1_iPSC <- e[[2]]

########## 2. patient WT HSC #########
# 2.1 Prepare mean_psi (nrows = sample) and expr_file (sample * gene name)
expr_psi_file = "/Users/chenxinyi/Desktop/output/data2_hHSC/HSC_p1_WT_expr_psi.csv"
expr_psi <- read.csv(expr_psi_file,header = TRUE,stringsAsFactors = FALSE, sep=" ") 
expr_file = SR_normed@assays$RNA@data
expr_file <- t(as.matrix(expr_file))
rownames(expr_file) <- gsub("TPM_","",rownames(expr_file))

#应该按照样本与SRRrun的对应关系组合一下，到学校再弄一下
