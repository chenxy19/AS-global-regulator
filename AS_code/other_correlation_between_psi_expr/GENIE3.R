"""
GENIE3: finding factors that best explain the response variable with Random Forest.

1. (old version) Run GENIE3 and rank all the X according to their relative importance.
   Input: expr_psi.csv from the lab cluster.
   Output: GENIE3_rank.csv
   But the ensmebl ID->gene name conversion is not effective here. See on the Seurat_HSC.R on cluster.
   
2. GO on HSC data
3. study genes annotated with nucleic acid activity
4. find n_feature and mean_psi correlation.
"""


########## 1. Run GENIE3 (old version) ##########
#BiocManager::install("GENIE3")
library(GENIE3)
HSC_p1_WT<- read.table("/Users/chenxinyi/Desktop/output/HSC_p1_WT_expr_psi.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 

nonERCC_index<-grep("^cxy",colnames(HSC_p1_WT),invert=TRUE)
#invert 可以直接选出没有这一pattern的!!!

HSC_p1_WT<-HSC_p1_WT[,nonERCC_index]

set.seed(123)
weightMat<-GENIE3(as.matrix(t(HSC_p1_WT)),targets="psi_HSC",nCores=1,nTrees=10000)  #1000是default
weight_index<-sort(weightMat, decreasing = TRUE, index.return=TRUE)
HSC_p1WT_GENIE3_rank<- colnames(HSC_p1_WT)[weight_index$ix]


write.table (HSC_p1WT_GENIE3_rank, file = "/Users/chenxinyi/Desktop/output/GENIE3_rank/HSC_p1WT_GENIE3_rank.csv", row.names =TRUE, col.names =TRUE, quote =TRUE)

# 到底有多少基因是ENSG开头的呀。。。三分之一都是ENSG，这可太不好了。。。

####### 2. GO: single gene effect ##########

top_genes <- colnames(HSC_p1_WT)[weight_index$ix[1:1000]]  

require(org.Hs.eg.db)  #这个包里存有人的注释文件
library(topGO)   #画GO图用的
library(clusterProfiler)
library(Rgraphviz)
library(pathview)

#将symbolID转换成ENTREZID
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = top_genes,
                       keytype = "SYMBOL",
                       column = "ENTREZID")   
#有67个基因没有转换成功
DEG.entrez_id = na.omit(DEG.entrez_id)
#biological process 富集分析
#gene_mat2<-read.csv (file = "/Users/chenxinyi/Desktop/gene_mat2.csv",sep=" ",header = TRUE,row.names = 1,stringsAsFactors = FALSE)
ref.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = colnames(HSC_p1_WT),
                       keytype = "SYMBOL",
                       column = "ENTREZID")   
#有779个基因没有转换成功(因为之前从ensemble转成基因名就已经失败了)
ref.entrez_id = na.omit(ref.entrez_id)

erich.go.BP_ref = enrichGO(gene = DEG.entrez_id,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID",
                           ont = "BP",
                           pvalueCutoff = 0.5,
                           qvalueCutoff = 0.5,
                           universe = ref.entrez_id)  # universe: background genes

##分析完成后，作图
dotplot(erich.go.BP_ref)

#没有什么显著的GO
########## 3. 看一下和核酸酶活性相关的基因 ############
GO_posit_regulat_nuclease_activity<-c("PCNA","HSPA1A","DDX11","TCEA1","RPS3","AKT1","HMGB2","PRKCD","psi_HSC")
GO_posit_regulat_nuclease_activity %in% colnames(HSC_expr_psi)
posit_regulat_nuclease_activity<-HSC_expr_psi[,GO_posit_regulat_nuclease_activity] 

GO_nega_regulat_nuclease_activity<-c("NEIL1","TERF1","DFFA","TMBIM6","ABCE1","TERF2","psi_HSC")
GO_nega_regulat_nuclease_activity %in% colnames(HSC_expr_psi)
nega_regulat_nuclease_activity<-HSC_expr_psi[,GO_nega_regulat_nuclease_activity] 

layout(matrix(1:1, 1, 1,byrow=T))
plotfunc<-function(GO){
  comb<-nega_regulat_nuclease_activity
  a<-cor.test(comb[,GO],comb$psi_HSC)
  print(a)
  print(GO)
  plot(comb[,GO],comb$psi_HSC,main=paste("Pearson corr=",round(a$estimate,5), "\np-value=",a$p.value,sep=""),xlab=GO,ylab="mean_psi (by sum of junctions)")
  LM<-lm(comb$psi_HSC ~ comb[,GO])
  abline(LM)
}
lapply(GO_nega_regulat_nuclease_activity[1:6], plotfunc)

GO_cell_division <- read.csv("cell_division.txt",header = FALSE,stringsAsFactors = FALSE) 
GO_DNA_damage <- read.csv("DNA_damage.txt",header = FALSE,stringsAsFactors = FALSE) 

find_gene_high_correlation<-function(comb,rowname_mean_psi,
                                     cor_method="pearson",
                                     cor_threshold=0.5){
  vol_data<-data.frame(matrix(NA,(ncol(comb)-1),2))
  rownames(vol_data)<-colnames(comb)[1:(ncol(comb)-1)]
  colnames(vol_data)<-c("corr","pvalue")
  for (i in 1:(ncol(comb)-1)){
    a<-cor.test(comb[,i],comb[,rowname_mean_psi],method = cor_method)
    vol_data[i,1] <-a$estimate
    vol_data[i,2] <-a$p.value
  }
  volcano<-subset(vol_data,select = c("corr","pvalue"))
  n_Bonferroni_correction=ncol(comb)-1
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
  return(e)
}

correlated_gene_pack<-find_gene_high_correlation(HSC_expr_psi,"psi_HSC",cor_method="spearman",cor_threshold=0.5)
posi_gene<-correlated_gene_pack[[1]]  # 476
nega_gene<-correlated_gene_pack[[2]]  # 42
correlated_gene_pack[[3]]

########### 4. n_feature & psi correlation ##################
#global mean psi
nfeature <- as.data.frame(SR_HSC$nFeature_RNA)
rownames(nfeature) <- gsub("NumReads_", "", rownames(nfeature))
HSC_expr_psi_nfeature <- merge_by_rownames(HSC_expr_psi,nfeature)
colnames(HSC_expr_psi_nfeature)[ncol(HSC_expr_psi_nfeature)]<-"nfeature"

LM <- lm(HSC_expr_psi_nfeature$psi_HSC ~ HSC_expr_psi_nfeature$nfeature)
a <- cor.test(HSC_expr_psi_nfeature$psi_HSC , HSC_expr_psi_nfeature$nfeature)
plot(HSC_expr_psi_nfeature$nfeature,HSC_expr_psi_nfeature$psi_HSC,
     xlab="n_features",ylab="mean psi (by adding up all junction reads)",
     main=paste("Pearson corr=",round(a$estimate,5), "\np-value=",a$p.value,sep=""))
abline(LM)
hist(HSC_expr_psi_nfeature$nfeature)

#gene-average psi

se_psi= read.csv("HSC2_se_bi_psi.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
se_psi_mean<-se_psi
se_psi_mean$mean_psi<-0
#小于10个细胞表达的就不要了

for (i in 1:length(se_psi[,1])){
  se_psi_mean$mean_psi[i]<-mean(as.numeric(se_psi[i,which(!is.na(se_psi[i,]))]))   #na.rm=TRUE
}

all_merged <- merge_by_rownames(HSC_expr_psi_nfeature,se_psi_mean)

#这个是gene average的
all_merged <- all_merged[which(!is.na(all_merged$mean_psi)),]
plot(all_merged$nfeature,all_merged$mean_psi)
LM<-lm(all_merged$mean_psi ~ all_merged$nfeature)
abline(LM)

all_merged_de01 <- all_merged[which(all_merged$mean_psi!=0 & all_merged$mean_psi!=1),]
a<-cor.test(all_merged_de01$mean_psi,all_merged_de01$nfeature)
plot(all_merged_de01$nfeature,all_merged_de01$mean_psi,xlab="n_feature",ylab="mean psi (gene-average psi)", main=paste("Pearson corr=",round(a$estimate,5), "\np-value=",a$p.value,sep=""))
LM<-lm(all_merged_de01$mean_psi ~ all_merged_de01$nfeature)
abline(LM)

#其实这种问题做GSEA更好，否则top基因数的选取会对结果有较大的影响
#然后把之前在前面三类细胞中也做一下GSEA。这样的GSEA和原先的GSEA有所不同，因为permutation只能对行来做，不能对列来做，这样就是彻底的随机化了，无法保留基因与基因之间的相关性。

#如果我直接去做GSEA分析呢？这种情况下做GSEA是最合适的，因为你可以直接evaluate这个gene set 是靠前还是靠后分布


