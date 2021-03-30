"""
GO analysis to find pathways associated with mean psi.

1. Identifying global factors related to gene-averaged psi. Plot the distribution of mean psi.
   Input: se_bi_psi.csv
   Output: histograms.

2. Violin Plot
   
3. GO enrichment analysis

4. correlation analysis and permutation to produce a background distribution
   
"""

########## 1. Identifying global factors ##########
# only bimodal AS genes histogram (can be ignored)
iPSC_se_psi= read.csv("/Users/chenxinyi/Desktop/output/iPSC_se_bi_psi.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
rownames(MN_se_psi)       #295和315是pool，295 QC时没给去掉，说明QC时总counts数还是要看的
iPSC_se_psi<-iPSC_se_psi[!rownames(iPSC_se_psi) %in% c("SRR4047260", "SRR4047278","SRR4047295","SRR4047315","SRR4047390","SRR4047332","SRR4047331", "SRR4047323"),] #NPC 31
iPSC_se_psi_mean<-iPSC_se_psi
MN_se_psi_mean$mean_psi<-0
MN_se_psi_mean$partial_mean20<-0
MN_se_psi_mean$partial_mean25<-0

cell_number_each_AS<-apply(iPSC_se_psi,2,sum_non_na)   #1 是对行进行操作

hist(as.numeric(cell_number_each_AS),main="number of cells each AS is detected in",xlab = "number of cells" )

which_gene_exceed_20 <- which(as.numeric(cell_number_each_AS) > 20)  #top 233/192
MN_se_psi_partial20<-MN_se_psi[,which_gene_exceed_20]
which_gene_exceed_25 <- which(as.numeric(cell_number_each_AS) > 25)  #top 70/79
MN_se_psi_partial25<-MN_se_psi[,which_gene_exceed_25]

for (i in 1:length(MN_se_psi[,1])){
  MN_se_psi_mean$mean_psi[i]<-mean(as.numeric(MN_se_psi[i,which(!is.na(MN_se_psi[i,]))]))   #na.rm=TRUE
  MN_se_psi_mean$partial_mean20[i]<-mean(as.numeric(MN_se_psi_partial20[i,which(!is.na(MN_se_psi_partial20[i,])) ]))
  MN_se_psi_mean$partial_mean25[i]<-mean(as.numeric(MN_se_psi_partial25[i,which(!is.na(MN_se_psi_partial25[i,])) ]))
}

hist(iPSC_se_psi_mean$mean_psi,col=rgb(1,0,0,1/4),breaks = 10,xlim = range(0.45,0.7),xlab = "Mean psi for all AS (gene-average)",main = "Histogram of mean psi across single iPSCs")
hist(iPSC_se_psi_mean$partial_mean20,add=T, col=rgb(0,0,1,1/4),breaks = 10 )
hist(iPSC_se_psi_mean$partial_mean25,add=T, col=rgb(0,1,0,1/4),breaks = 10)
legend("topright", c("all", "Top233","Top70"),fill=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4),rgb(0,1,0,1/4)),bty="n")  #bty="n"是legend周围的边框不要

######### 2. violin plot ###########

iPSC_se_psi_mean<-iPSC_se_psi

for (i in 0:8){
  iPSC_se_psi_mean[,ncol(iPSC_se_psi_mean)+1]<-0
  colnames(iPSC_se_psi_mean)[ncol(iPSC_se_psi_mean)] <- paste("psi_for_AS_detected_in_",as.character(5*i),"-",as.character(5*i+5),"_cells",sep="")
  which_gene_exceed <- which(as.numeric(cell_number_each_AS) > 5*i & as.numeric(cell_number_each_AS) <= 5*i+5)  
  iPSC_se_psi_partial<-iPSC_se_psi[,which_gene_exceed]
  
  for (row in 1:nrow(iPSC_se_psi)){
    iPSC_se_psi_mean[row,ncol(iPSC_se_psi_mean)]<-mean(as.numeric(iPSC_se_psi_partial[row,which(!is.na(iPSC_se_psi_partial[row,])) ]))
  }
}
library(reshape2)
library(ggplot2)
try_melt<-melt(iPSC_se_psi_mean[,tail(colnames(iPSC_se_psi_mean),9)],varnames=tail(colnames(iPSC_se_psi_mean),9))
try_melt$variable<-as.factor(try_melt$variable)
p<-ggplot(try_melt, aes(x = variable, y = value)) #注释：”x=”，”y=”表示x轴和y轴表示的变量数值，p表示图像对象
p+geom_violin() #注释：画出violin plot的函数
p+geom_violin(aes(fill = variable))+  labs(title = "iPSC", x="Dose (mg)", y = "mean psi for all iPSC")+theme(plot.title = element_text(hjust = 0.5))  #标题居中

plot(8*(0:3),apply(iPSC_se_psi_mean[,tail(colnames(iPSC_se_psi_mean),4)],2,mean),xlab="number of cells in which the AS is detected",ylab="mean psi in the bin (gene-average)")

#all genes histogram
NPC_se_psi_all= read.csv("/Users/chenxinyi/Desktop/NPC_psi/se/psi.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
rownames(MN_se_psi_all)       #295和315是pool，295 QC时没给去掉，说明QC时总counts数还是要看的
filter116= read.csv("/Users/chenxinyi/Desktop/filtered116",sep=" ",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 

MN_se_psi_all2<-MN_se_psi_all2[which(rownames(MN_se_psi_all) %in% filter116$x),]  #MN 61->48, NPC 32, iPSC 33
#NPC_se_psi_all<-NPC_se_psi_all[!rownames(NPC_se_psi_all) %in% c("SRR4047260","SRR4047278","SRR4047390"),] #NPC 31
#MN_se_psi_all<-MN_se_psi_all[!rownames(MN_se_psi_all) %in% c("SRR4047295","SRR4047315"),] #NPC 31
MN_se_psi_mean2<-MN_se_psi_all2
MN_se_psi_mean2$mean_psi<-0
for (i in 1:length(MN_se_psi_all2[,1])){
  MN_se_psi_mean2$mean_psi[i]<-mean(as.numeric(MN_se_psi_all2[i,which(!is.na(MN_se_psi_all2[i,]))]))   #na.rm=TRUE
}
hist(MN_se_psi_mean$mean_psi,breaks = 10,xlab = "mean_psi",main = "Histogram of mean psi across single MNs")

#决定用后者（All genes 得到的分类结果来做differential analysis）
SR<-readRDS("/Users/chenxinyi/Desktop/filtered116.rds", refhook = NULL)
SR@meta.data$sample_id<-gsub("TPM_", "", rownames(SR@meta.data))

low<-rownames(NPC_se_psi_mean)[which(NPC_se_psi_mean$mean_psi<0.62)]
high<-rownames(NPC_se_psi_mean)[which(NPC_se_psi_mean$mean_psi>=0.63)]

#要初始化标签！不然会跟之前的标签搞混！！！
Idents(object=SR)<-0
Idents(object = SR, cells = which(SR@meta.data$sample_id %in% low)) <- 'NPC_low'
Idents(object = SR, cells = which(SR@meta.data$sample_id %in% high)) <- 'NPC_high'

low_markers <- FindMarkers(SR, ident.1 = "MN_low", ident.2 = "MN_high", min.pct = 0.25,only.pos=TRUE)   
low_markers_NPC <- FindMarkers(SR, ident.1 = "NPC_low", ident.2 = "NPC_high", min.pct = 0.25,only.pos=TRUE)   
low_markers_iPSC <- FindMarkers(SR, ident.1 = "iPSC_low", ident.2 = "iPSC_high", min.pct = 0.25,only.pos=TRUE)   
#low psi markers MN 822, NPC 661, iPSC 100
high_markers <- FindMarkers(SR, ident.1 = "MN_high", ident.2 = "MN_low", min.pct = 0.25,only.pos=TRUE)   
high_markers <- FindMarkers(SR, ident.1 = "NPC_high", ident.2 = "NPC_low", min.pct = 0.25,only.pos=TRUE)   
high_markers_iPSC <- FindMarkers(SR, ident.1 = "iPSC_high", ident.2 = "iPSC_low", min.pct = 0.25,only.pos=TRUE)   
# high markers MN 954, NPC 832, iPSC 166
# min.pct: speed up by only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations.
# 有警告：wilcox.test.default: cannot compute exact p-value with ties
# only.pos=FALSE by default，就是多的少的都会有
low_markers<-intersect(rownames(low_markers_MN),rownames(low_markers_NPC))  #71
low_markers<-intersect(a,rownames(low_markers_iPSC)) #三个的交集只有1个了
#MN, NPC 与iPSC相交的只有6个和7个

high_markers<-intersect(rownames(high_markers_MN),rownames(high_markers_NPC))  #135
high_markers<-intersect(high_markers,rownames(high_markers_iPSC)) #only 6

head(low_markers, n = 6)

VlnPlot(SR, idents=c("MN_low","MN_high"),features = rownames(low_markers)[1:6])+scale_color_discrete()

######### 3. GO enrichment analysis ############

#BiocManager::install("org.Hs.eg.db")
require(org.Hs.eg.db)  #这个包里存有人的注释文件
#BiocManager::install("topGO")
library(topGO)   #画GO图用的
#BiocManager::install("clusterProfiler")  #用来做富集分析
#BiocManager::install("Rgraphviz")
#BiocManager::install("pathview") #看KEGG pathway的
library(clusterProfiler)
library(Rgraphviz)
library(pathview)

 #获得基因 symbol ID
#将symbolID转换成ENTREZID
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = shared_weight_ordered,
                       keytype = "SYMBOL",
                       column = "ENTREZID")   
#有67个基因没有转换成功
DEG.entrez_id = na.omit(DEG.entrez_id)
#biological process 富集分析
#gene_mat2<-read.csv (file = "/Users/chenxinyi/Desktop/gene_mat2.csv",sep=" ",header = TRUE,row.names = 1,stringsAsFactors = FALSE)
ref.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = rownames(SR@assays$RNA@data),
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

erich.go.CC_ref = enrichGO(gene = DEG.entrez_id,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)#,
                       #universe = ref.entrez_id)
## 画图
barplot(erich.go.CC_ref)

erich.go.MF_ref = enrichGO(gene = DEG.entrez_id,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID",
                           ont = "MF",
                           pvalueCutoff = 0.5,
                           qvalueCutoff = 0.5)#,
                           #universe = ref.entrez_id)  # universe: background genes
barplot(erich.go.MF_ref)

erich.KEGG_ref=enrichKEGG(gene=DEG.entrez_id,
                          organism = "hsa",
                          keyType = "kegg",
                          pvalueCutoff = 0.5,
                          qvalueCutoff = 0.5,
                          universe = ref.entrez_id)
barplot(erich.KEGG_ref)


############# 4. 求correlation ###############
normed_expr<-as.data.frame(SR_normed@assays$RNA@data) #将一个稀疏矩阵转化成正常的矩阵
normed_expr<-t(normed_expr)
rownames(normed_expr)<-gsub("TPM_", "", rownames(normed_expr))

iPSC_expr_psi <- merge(normed_expr,iPSC_global_mean_psi,by=0,all=FALSE )
#by=0和“row.names”都是用行名来做合并
iPSC_expr_psi<-iPSC_expr_psi[,-1]
iPSC_expr_psi2<-iPSC_expr_psi[,-((length(colnames(iPSC_expr_psi))-3):(length(colnames(iPSC_expr_psi))-1))]
#MN_expr_psi2[,grep("^isoform", colnames(MN_expr_psi2))]  #R正则表达式,检查一下isoform删干净了
#第13676列是mean_psi
#删掉全是0的基因；前面normalization是在所有细胞上做的，这里只有一种细胞
iPSC_expr_psi2<-iPSC_expr_psi2[,which(colSums(iPSC_expr_psi2)!=0)]

row<- data.frame(matrix(NA,2,length(colnames(iPSC_expr_psi2))-1))
rownames(row)<-c("pearson_cor","spearman_cor")
colnames(row)<-colnames(iPSC_expr_psi2)[1:(length(colnames(iPSC_expr_psi2))-1)]

for (i in 1:(length(colnames(iPSC_expr_psi2))-1)){
  row[1,i] <-cor(iPSC_expr_psi2[,i],iPSC_expr_psi2[,length(colnames(iPSC_expr_psi2))],method = "pearson")
  row[2,i] <-cor(iPSC_expr_psi2[,i],iPSC_expr_psi2[,length(colnames(iPSC_expr_psi2))],method = "spearman")
  }
hist(as.numeric(row[1,]),breaks = 10,xlab = "pearson corr",main = "Histogram of pearson corr between mean psi \nand gene expression across single NPCs")

########### 4.2 Permutation ###############
iter<-300
row2<- data.frame(matrix(NA,iter,length(colnames(iPSC_expr_psi2))-1))
#rownames(row2)<-"pearson_cor"
colnames(row2)<-colnames(iPSC_expr_psi2)[1:(length(colnames(iPSC_expr_psi2))-1)]

for (j in 1:iter){
  rand<-sample(iPSC_expr_psi2[,length(colnames(iPSC_expr_psi2))])
  for (i in 1:(length(colnames(iPSC_expr_psi2))-1)){
    row2[j,i] <-cor(iPSC_expr_psi2[,i],rand,method = "spearman")
  }
}
library(reshape2)
agg<-melt(row2)
hist(as.numeric(row[2,]),freq=FALSE,col=rgb(1,0,0,1/4),breaks = 20,xlim=range(-1,1),xlab = "spearman corr",add=T,main = "Histogram of spearman corr between randomized psi \nand gene expression across single iPSCs")
hist(agg$value,col=rgb(0,0,1,1/4),freq=FALSE,breaks = 40,xlim=range(-1,1),xlab = "spearman corr",main = "Histogram of spearman corr between mean psi \nand gene expression across single iPSCs")
legend("topright", c("iPSC", "permutated"),fill=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)),bty="n")  #bty="n"是legend周围的边框不要

t.test(as.numeric(row[2,]),agg$value)



posi_cor<-row[which(row[2,]>0.4)]
nega_cor<-row[which(abs(row[2,]) > 0.4 & row[1,]<0)]
posi_top<-colnames(posi_cor)[order(posi_cor[1,],decreasing=TRUE)]  #进行排序
nega_top<-colnames(nega_cor)[order(nega_cor[1,],decreasing=FALSE)]  #进行排序

iPSC_posi_top<-posi_top
iPSC_nega_top<-nega_top

NPC_posi_top<-posi_top
NPC_nega_top<-nega_top

#row[1,order(row)[13500:13520]]    #order可以对dataframe 进行排序
MN_posi_top<-posi_top
MN_nega_top<-nega_top

