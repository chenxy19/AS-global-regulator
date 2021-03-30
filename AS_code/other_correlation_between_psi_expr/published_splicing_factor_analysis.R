"""
Analyze whether well studied splicing factors explain the variation in mean psi.
Input: all_mean_psi_corrected0831.csv

1. Correlation and GENIE3 analysis with mean psi.
2. 
"""

expr_mat<-as.matrix(SR_normed@assays$RNA@data)
splic_factor <- c("SRSF1","SRSF2","SRSF3","SRSF4","SRSF5","SRSF6","SRSF7","SRSF8","SRSF9",
  "SRSF10","TRA2A", "SRRM1", "HNRNPA1", 
  "HNRNPA2B1","HNRNPC", "HNRNPF", "HNRNPH1", "HNRNPL","PTBP1", "PTBP2", "SFPQ", "RBM4", "RBM5",
  "RBFOX2", "RBM14", "RBM17", "RBM23", "RBM25", "RBM39","RBFOX1", "MBNL1", "MBNL2", 
  "MBNL3", "CELF1", "CELF2", "CELF4", "QKI", "KHDRBS3","KHDRBS1", "NOVA1"," NOVA2",
  "ELAVL4", "TIA1", "TIAL1") 
# 为啥表达矩阵里会没有这些基因。。。原来这些都是别称。。。
                          
SF_expr <- expr_mat[which(rownames(expr_mat) %in% splic_factor),]             
SF_expr<-as.data.frame(t(SF_expr))                          
rownames(SF_expr)<-gsub("TPM_", "", rownames(SF_expr))

all_psi_new= read.csv("/Users/chenxinyi/Desktop/output/all_mean_psi_corrected0831.csv",sep=" ",header = TRUE,stringsAsFactors = FALSE) 
#坏处是这里只有111个基因，有时间可以重新把别的细胞也求一下
all_SF_psi<-merge_by_rownames(SF_expr,all_psi_new)

#先画一下热图
#BiocManager::install("ComplexHeatmap")
library(gplots)
library(ComplexHeatmap)
col <- colorRampPalette(c("red", "white", "blue"))(256)
Heatmap(as.matrix(all_SF_psi[,-(43:45)]),show_row_names = FALSE,name = "expression", col=bluered(100), split = all_SF_psi$cell_type )

########## 1. Correlation analysis with global mean psi ##########
#find_gene_high_correlation()在recalculate_psi_less_coverage.R
correlated_gene_pack<-find_gene_high_correlation(all_SF_psi,"glbl_mean_psi",cor_method="spearman",cor_threshold=0.5)
posi_gene<-correlated_gene_pack[[1]]
nega_gene<-correlated_gene_pack[[2]]
correlated_gene_pack[[3]]  #作图

library(ggpmisc)

SF_list<-c("ELAVL4","NOVA1","RBFOX1","CELF4","PTBP1","SRSF7","SRSF1","HNRNPA","HNRNPF","HNRNPA2B1")
gene_mean_psi_scatter(SF_list[5])
#list_of_plot<-lapply(SF_list,gene_mean_psi_scatter)

######### 2. GENIE3 find regulatory relationship with RF ##########
#应该只找在MN中highly variable的基因！
#不挑highly variable gene看一看？

Idents(object = SR, cells = 1:62) <- 'MN'
Idents(object = SR, cells = 63:98) <- 'NPC'
Idents(object = SR, cells = 99:141) <- 'iPSC'
SR[["percent.mt"]] <- PercentageFeatureSet(SR, pattern = "^MT-")
#SR_iPSC<-subset(SR, idents="NPC")

#SR_qc <- subset(SR, subset = nFeature_RNA < 12000 & nFeature_RNA > 2000 & percent.mt < 15)
SR_qc<-SR[,!colnames(SR) %in% c("TPM_SRR4047260", "TPM_SRR4047278","TPM_SRR4047295","TPM_SRR4047315","TPM_SRR4047390","TPM_SRR4047332","TPM_SRR4047331", "TPM_SRR4047323")] 
SR_normed<- NormalizeData(SR_qc, normalization.method = "LogNormalize", scale.factor = 10000)
#SR_normed <- FindVariableFeatures(SR_normed, selection.method = "vst", nfeatures = 2000)

expr_mat<-as.matrix(SR_normed@assays$RNA@data)
colnames(expr_mat)<-gsub("TPM_", "", colnames(expr_mat))

#expr_mat_var<-expr_mat[SR_normed@assays$RNA@var.features,]  #5000个variable genes
expr_mat_var<-t(expr_mat)
#一共有13674个基因

library(GENIE3)
set.seed(123)
expr_psi<-merge_by_rownames(expr_mat_var,all_psi_new) #最后三列是all_psi_new里面的

single_type_expr_psi <- t(expr_psi[,-((ncol(expr_psi)-1):ncol(expr_psi))])
weightMat<-GENIE3(as.matrix(single_type_expr_psi),targets="glbl_mean_psi",nCores=1,nTrees=10000)  #1000是default

weightmat<-as.data.frame(weightMat)
weight_index<-sort(weightmat$glbl_mean_psi, decreasing = TRUE, index.return=TRUE)
weight_ordered_gene<-rownames(weightmat)[weight_index$ix]

iPSC_weight_ordered<-weight_ordered_gene
MN_weight_ordered<-weight_ordered_gene
NPC_weight_ordered<-weight_ordered_gene
top_number<-1000
a<-intersect(iPSC_weight_ordered[1:top_number],MN_weight_ordered[1:top_number])
shared_weight_ordered<-intersect(a,NPC_weight_ordered[1:top_number])

# top 1000有100个intersect的
GO_posit_regulat_nuclease_activity<-c("PCNA","HSPA1A","DDX11","TCEA1","RPS3","AKT1","HMGB2","PRKCD")
GO_posit_regulat_nuclease_activity %in% rownames(SR@assays$RNA@data)
#都在！太好了！
#只保留8列
posit_regulat_nuclease_activity<-expr_mat_var[,GO_posit_regulat_nuclease_activity] 
expr_psi<-merge_by_rownames(posit_regulat_nuclease_activity,all_psi_new) #最后三列是all_psi_new里面的

GO_nega_regulat_nuclease_activity<-c("NEIL1","TERF1","DFFA","TMBIM6","ABCE1","TERF2")
GO_nega_regulat_nuclease_activity %in% rownames(SR@assays$RNA@data)
nega_regulat_nuclease_activity<-expr_mat_var[,GO_nega_regulat_nuclease_activity] 
expr_psi<-merge_by_rownames(nega_regulat_nuclease_activity,all_psi_new) #最后三列是all_psi_new里面的

GO_RPS<-colnames(expr_psi)[grep("^RPS",colnames(expr_psi))]
RPS<-expr_mat_var[,GO_RPS] 
expr_psi<-merge_by_rownames(RPS,all_psi_new) #最后三列是all_psi_new里面的



col <- colorRampPalette(c("red", "white", "blue"))(256)
#Heatmap(as.matrix(expr_psi[,1:8]),show_row_names = FALSE,name = "expression", col=bluered(100), split = expr_psi$cell_type )
sort_index<-sort(expr_psi$RPS3, index.return = TRUE)
Heatmap(as.matrix(expr_psi[,1:8]),show_row_names = FALSE,name = "expression", col=bluered(100), row_order=sort_index$ix)
row_order=expr_psi$glbl_mean_psi 
# row_order = NULL,
#column_order = NULL,

plotfunc<-function(comb,GO,x=GO){
  a<-cor.test(comb[,GO],comb$glbl_mean_psi)
  print(a)
  plot(comb[,GO],comb$glbl_mean_psi,xlab=x,ylab="mean_psi (by sum of junctions)",main=paste("Pearson corr=",round(a$estimate,5), "\np-value=",a$p.value,sep=""))
  LM<-lm(comb$glbl_mean_psi ~ comb[,GO])
  abline(LM)
}
a<-rep(0,44)
for (i in 1:44){
  a[i]<-cor(expr_psi[,i],expr_psi$glbl_mean_psi)
}
which(abs(a)>0.5)
RPS_name<-colnames(expr_psi)[which(abs(a)>0.5)][7]
RPS_name

ggplot(as.data.frame(expr_psi),aes(x=`RPS15A`,y=`glbl_mean_psi`,color=as.factor(cell_type))) + geom_point(shape=19) +
  geom_point() + 
  scale_colour_brewer(palette = "Set1") +
  geom_smooth(method="lm",show.legend = TRUE,level=0.95) +
  #scale_fill_discrete(breaks=c("5_n1","5_n2","3_n1","3_n2")) +
  xlab("RPS20 expression level") + ylab("mean psi (by summing up junction reads)") 

layout(matrix(1:4, 2, 2,byrow=T))
plotfunc(expr_psi[which(expr_psi$cell_type=="MN"),],RPS_name)
plotfunc(expr_psi[which(expr_psi$cell_type=="NPC"),],RPS_name)
plotfunc(expr_psi[which(expr_psi$cell_type=="iPSC"),],RPS_name)
plotfunc(expr_psi[,],RPS_name)

cell_division_gene= read.csv("/Users/chenxinyi/Desktop/cell_division.txt",header = FALSE,stringsAsFactors = FALSE) 
cell_division_gene<-as.vector(cell_division_gene$V1)

high_cor_cell_division<- c("CTDP1","TUBG1","FIGNL1","DYNLT1","RPS3", "CDC123","WNT9B","CDK2AP2","BAG6","RACK1")
cell_division<-expr_mat_var[,which(colnames(expr_mat_var) %in% high_cor_cell_division)] 
expr_psi<-merge_by_rownames(cell_division,all_psi_new) #最后三列是all_psi_new里面的
expr_psi_iPSC<-expr_psi[which(expr_psi$cell_type=="iPSC"),]

layout(matrix(1:8, 2, 4,byrow=T))
for (i in 9:16){
  plotfunc(expr_psi_iPSC,colnames(expr_psi_iPSC)[i])
}




find_gene_high_correlation<-function(comb,rowname_mean_psi,
                                     cor_method="pearson",
                                     cor_threshold=0.5){
  vol_data<-data.frame(matrix(NA,(ncol(comb)-3),2))
  rownames(vol_data)<-colnames(comb)[1:(ncol(comb)-3)]
  colnames(vol_data)<-c("corr","pvalue")
  for (i in 1:(ncol(comb)-3)){
    a<-cor.test(comb[,i],comb[,rowname_mean_psi],method = cor_method)
    vol_data[i,1] <-a$estimate
    vol_data[i,2] <-a$p.value
  }
  volcano<-subset(vol_data,select = c("corr","pvalue"))
  n_Bonferroni_correction=ncol(comb)-3
  pvalue_threshold<-0.05/n_Bonferroni_correction
  threshold<-as.factor((abs(volcano$corr)>cor_threshold) & volcano$pvalue<pvalue_threshold)   #其实最好是用bonferonni corrected P value
  r03=ggplot(volcano,aes(corr,-log2(pvalue),colour=threshold))+geom_point()
  r04=r03+labs(title=paste("Volcanoplot of ",cor_method," corr with mean psi\n(Bonferroni-corrected test)",sep=""))+theme(plot.title = element_text(hjust = cor_threshold))+xlim(-1,1)
  r05=r04+geom_vline(xintercept=c(-1*cor_threshold,cor_threshold),linetype="dotted",size=1)+geom_hline(yintercept=-log2(pvalue_threshold),col="blue")
  
  library(ggrepel)    #可以让字与字自动错开
  label<-rownames(volcano)
  label[which((abs(volcano$corr)<=cor_threshold) | volcano$pvalue>=pvalue_threshold)]<-NA
  r06=r05+geom_text_repel(aes(corr, -log2(pvalue), label = label))
  r06
  
  posi_gene<-label[which(volcano$corr>cor_threshold & volcano$pvalue<pvalue_threshold)]
  nega_gene<-label[which(abs(volcano$corr)>cor_threshold & volcano$corr<0 & volcano$pvalue<pvalue_threshold)]
  e <- list(posi_gene,nega_gene,r06)  #因为函数只能返回一个值，所以给打包成一个对象
  return(vol_data)
}





expr_psi$positive_activity <- apply(expr_psi[,1:8],1,sum)
  
plotfunc(expr_psi[which(expr_psi$cell_type=="iPSC"),],"positive_activity","positive activity in iPSC")
plotfunc(expr_psi[which(expr_psi$cell_type=="MN"),],"positive_activity","positive activity in MN")
plotfunc(expr_psi[which(expr_psi$cell_type=="NPC"),],"positive_activity","positive activity in NPC")

# rapid!
# 
#weight_ordered_gene, SR_normed@assays$RNA@var.features