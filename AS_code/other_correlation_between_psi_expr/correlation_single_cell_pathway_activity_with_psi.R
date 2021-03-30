"""
Analyze whether the activity of single pathway activity in single cells are correlated to mean psi.

1. Genarate regulon
   Output: dimred/GO-BP-Hs-MSigDB-regulon.rda

2. prepare gene expr matrix and load psi
   Input: Seurat object SR_normed from Seurat.R, se_bi_psi.csv
   Output: expr_mat

3. Permutation to generate a backgorund distribution for correlation between expression and mean psi.

4. Set threshold and find shared significant gene sets across cell types.
"""

BiocManager::install("viper")
install.packages("msigdbr")

############ 1. Generate regulon ##############
#https://github.com/hd2326/BiologicalProcessActivity/blob/master/MSigDB-regulon/MSigDB-regulon.R#L7
library(viper)
library(msigdbr)
#generate GO-BP-Hs-MSigDB regulon

table <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")#GO, BP
#A data frame of gene sets with one gene per row.(662,851 x 9)
#category和subcategory的限制是可以不加的
# msigdbr_show_species()人和鼠都可以用

len <- table(table$gs_name)
#显示每个geneset有多少基因,长度7530，说明有这么多的geneset

gset <- lapply(names(len), function(x, table){
  gene <- table$gene_symbol[grep(x, table$gs_name)]
  list(tfmode=structure(rep(1, length(gene)), names=gene), likelihood=rep(1, length(gene)))
}, table=table)
#lapply是将第二个函数应用在第一个变量上；grep(x, table$gs_name)是找出每个gene set对应基因在table中的序号,从而提取出基因名,注意lapply将会返回一个list
#take a very lomng time，形成一个data structure

names(gset) <- names(len)
save(gset, file = "/Users/chenxinyi/Desktop/dimred/GO-BP-Hs-MSigDB-regulon.rda")

load(file="/Users/chenxinyi/Desktop/output/dimred/GO-BP-Hs-MSigDB-regulon.rda")
#不需要把load的结果附给gset，直接就可以了

########### 2. prepare gene expr matrix  ############
#直接读那个做好了normalization和对数转换的 SR_normed 就可以
expr_mat<-SR_normed@assays$RNA@data  #如果把稀疏矩阵展开成普通矩阵的话会非常非常占内存，一般的计算用这个线性的数据结构也能进行计算，只是慢一些而已
#normalize后的数据,normalize确保每个细胞的测序深度是一样的，是不能省的步骤

#注意：用循环超级慢的，用apply会非常快，但是回复的数据类型受到限制

len <- unlist(lapply(gset, function(x) length(x$tfmode)))  #和前面的len意思是一样的
gset <- gset[len >= 5 & len <= 500]   #根据genes set 的大小来过滤一下gene set,原先是50-100，试试放宽限制
nes_Hs <- aREA(as.matrix(expr_mat), gset)$nes  #as.matrix 可以将稀疏矩阵化作正常矩阵

nes_Hs<-t(nes_Hs)
rownames(nes_Hs)<-gsub("TPM_", "", rownames(nes_Hs))

######## find highly variable pathway activity and differential activity across #############
sc_activ <- CreateSeuratObject(counts = t(nes_Hs), project = "sc_activ_project")
#124个sample,3139个pathway activity
#标记一下cell type
Idents(object = sc_activ, cells=1:55) <- 'MN'
Idents(object = sc_activ, cells=56:90) <- 'NPC'
Idents(object = sc_activ, cells=91:124) <- 'iPSC'

sc_activ <- FindVariableFeatures( sc_activ, selection.method = "mvp", nfeatures = 2000)
# vst不足以用一个多项式关系normalize掉variance与表达量的关系，所以用mvp，用分bin的方法对每个bin分别求z score来标准化，坏处是没法控制nfeature
# mvp: The purpose of this is to identify variable features while controlling for the strong relationship between variability and average expression.
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sc_activ ), 10)  #选择variable feature的方法对后面影响很大

# plot variable features with and without labels
VariableFeaturePlot(sc_activ)
LabelPoints(plot = plot1, points = top10,repel = TRUE,xnudge =0, ynudge=0)  #repel是为了字和字之间的间隔

MN_markers <- FindMarkers(sc_activ, ident.1 = "MN", min.pct = 0.25,only.pos=FALSE)   
NPC_markers <- FindMarkers(sc_activ, ident.1 = "NPC", min.pct = 0.25,only.pos=FALSE)   
iPSC_markers <- FindMarkers(sc_activ, ident.1 = "iPSC", min.pct = 0.25,only.pos=FALSE)   
a<-union(rownames(MN_markers),rownames(NPC_markers))
Dpath<-union(a,rownames(iPSC_markers))   #一共449 个DEG

#gene_pool<-setdiff(rownames(sc_activ@assays$RNA@counts),Dpath)   #如果一定要在highly variable gene中选的话，就剩3个了。。。
#这样剩下2690

#先不过滤试下


############ load psi ################
#已经filtered过的细胞和基因
NPC_global_mean_psi= read.csv("/Users/chenxinyi/Desktop/output/NPC_global_mean_psi.csv",sep=" ",header = TRUE,stringsAsFactors = FALSE) 

NPC_se_psi= read.csv("/Users/chenxinyi/Desktop/output/NPC_se_bi_psi.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
NPC_se_psi_mean<-NPC_se_psi
NPC_se_psi_mean$mean_psi<-0
for (i in 1:length(NPC_se_psi[,1])){
  NPC_se_psi_mean$mean_psi[i]<-mean(as.numeric(NPC_se_psi[i,which(!is.na(NPC_se_psi[i,]))]))   #na.rm=TRUE
}
hist(NPC_se_psi_mean$mean_psi,breaks = 10,xlab = "mean_psi",main = "Histogram of mean psi across single iPSCs")

hist(MN_global_mean_psi$mean_psi_glbl,breaks = 20,xlab = "mean_psi",main = "Histogram of mean psi of MNs (new)")

NPC_DR<-merge(nes_Hs,NPC_global_mean_psi,by=0,all=FALSE) #48行是细胞数，4299 是isoform和gene set数
#先根据行名来merge

NPC_DR<-NPC_DR[,-1]
NPC_DR<-NPC_DR[,-((length(NPC_DR[1,])-3):(length(NPC_DR[1,])-1))]

row<- data.frame(matrix(NA,1,848))
rownames(row)<-"spearman_cor"
colnames(row)<-colnames(NPC_DR)[1:848]

for (i in 1:848){
  row[1,i] <-cor(NPC_DR[,i],NPC_DR[,849],method = "spearman")
}

################ 3. Permutation ##################
#set.seed(42)
iter<-1000
row2<- data.frame(matrix(NA,iter,848))
#rownames(row2)<-"pearson_cor"
colnames(row2)<-colnames(NPC_DR)[1:848]

for (j in 1:iter){
  rand<-sample(NPC_DR[,849])
  for (i in 1:848){
    row2[j,i] <-cor(NPC_DR[,i],rand,method = "spearman")
  }
}
library(reshape2)
agg<-melt(row2)
hist(as.numeric(row[1,]),freq=FALSE,col=rgb(1,0,0,1/4),breaks = 20,xlim=range(-1,1),xlab = "spearman corr",main = "Histogram of spearman corr between randomized psi \nand pathway activity across single NPCs")
hist(agg$value,col=rgb(0,0,1,1/4),freq=FALSE,breaks = 40,xlim=range(-1,1),xlab = "spearman corr",add=T,main = "Histogram of spearman corr between mean psi \nand pathway activity across single NPCs")
legend("topright", c("NPC", "permutated"),fill=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)),bty="n")  #bty="n"是legend周围的边框不要

t.test(as.numeric(row[1,]),agg$value)

######### 4. Test for each gene and find gene sets for regression ##########
#要求1）对至少两类细胞 2）significance 在95% 3）corr的绝对值大于0.4
row3<-t(row2)
colnames(row3)<-c(1:iter)
row3<-as.data.frame(row3)
row3$signif <-0
for (i in 1:length(row3[,1])){
#  sht<-shapiro.test(row2[,i])   #p>0.05，符合正态分布;0.074近似正态分布
#  if (sht$p.value<=0.05){
#    print(rownames(row3)[i])
#    print("not a normal distribution")   #大部分都不正态
#  }
  signif<-0.95
  low<-quantile(row2[,i],(1-signif)/2)
  up<-quantile(row2[,i],(1+signif)/2)
  if ( row[,i]>up | row[,i]<low){
    row3$signif[i]<-1
  }
}

posi_cor<-colnames(row)[which(row[1,]>0.4)]
nega_cor<-colnames(row)[which(abs(row[1,]) > 0.4 & row[1,]<0)]
posi_cor<-t(row[posi_cor])
posi_top<-rownames(posi_cor)[order(posi_cor,decreasing=TRUE)]  #进行排序
nega_cor<-t(row[nega_cor])
nega_top<-rownames(nega_cor)[order(nega_cor,decreasing=FALSE)]#[1:10]


sig_list<-rownames(row3)[which(row3$signif==1)]
iPSC_posi_list<-intersect(sig_list,posi_top) #4个
iPSC_nega_list<-intersect(sig_list,nega_top) #9个

NPC_posi_list<-intersect(sig_list,posi_top) #24个
NPC_nega_list<-intersect(sig_list,nega_top) #44个

MN_posi_list<-intersect(sig_list,posi_top) #35个
MN_nega_list<-intersect(sig_list,nega_top) #69个

a<-intersect(iPSC_posi_list,NPC_posi_list)
posi_all<-intersect(a,MN_posi_list)

b<-intersect(iPSC_nega_list,NPC_nega_list)
nega_all<-intersect(b,MN_nega_list)
