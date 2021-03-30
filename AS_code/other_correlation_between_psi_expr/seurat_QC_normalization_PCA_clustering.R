"""
Single cell data into a Seurat object.

Seurat single cell data analysis
   Input: all_gene_TPM.csv
   1. QC
   2. Normalization
   3. highly variable genes
   4. PCA
   5. clustering
   6. correlation analysis, but have better version in recalculate_psi_less_coverage.R
"""

######## *2. the second step in aggregating sailfish result, can be skipped now ############
cell_type="iPSC"
library(dplyr)
raw_folder<-paste("/Users/chenxinyi/Desktop/",cell_type,"_filted_fa",sep="")
SRR_list<-list.files(raw_folder)
for (folder in SRR_list){
  SRR_file=paste(raw_folder,folder,"quant.sf",sep="/")
  tpm= read.csv(SRR_file,header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
  gene_tpm <- group_by(tpm, ENSG) %>% summarize(TPM = sum(TPM))
  out_file=paste(raw_folder,"/",folder,"/",cell_type,"_gene_TPM.txt",sep="")
  write.table (gene_tpm, file = out_file, row.names =TRUE, col.names =TRUE, quote =TRUE)
}
#Groups:   ENSG [20,738]. group_by will not change how the data is arranged
#行名自动变成ENSG，TPM求和

######## Single cell analysis ############
# 都合并到一起分析
all_mat= read.csv("/Users/chenxinyi/Desktop/output/all_gene_TPM.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
#20738*141
#把ensemble ID (和GENCODE编号很像)和基因名称对应起来（因为需要通过MT来判断mitochondria gene）
library("biomaRt")
enID<-rownames(all_mat)
enID<-as.character(enID)
enID <- sub("[.][0-9]*","",enID)
all_mat$gene_id<-enID
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = enID, mart= mart)
gene_mat<-left_join(all_mat, genes, by = c("gene_id"="ensembl_gene_id"))  #类似于merge,左边的gene_id对应着右边的gene_name
#处理掉重复项,NA 和没有对上号的，用ENSG代替基因名
gene_mat[which(gene_mat$hgnc_symbol==""),"hgnc_symbol"]<-gene_mat[which(gene_mat$hgnc_symbol==""),"gene_id"]
gene_mat[which(is.na(gene_mat$hgnc_symbol)),"hgnc_symbol"]<-gene_mat[which(is.na(gene_mat$hgnc_symbol)),"gene_id"]

gene_mat[19967,"hgnc_symbol"]<-"PINX1.1"
#gene_mat[which(gene_mat$hgnc_symbol=="PINX1"),]  #19674 19968
#x<-gene_mat$hgnc_symbol
#y<-x[!duplicated(x)]    #感觉重复的超级多。。。。

rownames(gene_mat)<-gene_mat$hgnc_symbol

############ QC of genes ################
#genes with TPM>1 in at least 10 cells are identified according to the paper
gene_mat<-gene_mat[,-(142:143)]
#20739*141
gene_mat2<-gene_mat
#for (i in length(gene_mat[,1]):1){
#  if (sum(gene_mat2[i,]>1) < 10){     #按照文章把低表达量的基因给滤掉了，有可能严苛了，可以调整的
#    print(i)
#    gene_mat2<-gene_mat2[-i,]
#  }
#}
# 20740->13674；因为在这一步一些低表达量基因被剔除了，所以TPM总和不再是整的
write.table (gene_mat2, file = "/Users/chenxinyi/Desktop/gene_mat2.csv")

gene_mat2= read.csv("/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/gene_mat2.csv",sep=" ",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
library(Seurat)
SR <- CreateSeuratObject(counts = gene_mat2, project = "neuronAS_project")
SR
#An object of class Seurat 
#13674 features across 141 samples within 1 assay 
#Active assay: RNA (13674 features)

#加入cell type的label,和放进来之前的顺序是一致的 (往meta.data里面加，就是一个正常的dataframe)
Idents(object = SR, cells = 1:62) <- 'MN'
Idents(object = SR, cells = 63:98) <- 'NPC'
Idents(object = SR, cells = 99:141) <- 'iPSC'

############ 1. QC of cells ################
#the set of all genes starting with MT- as a set of mitochondrial genes
SR[["percent.mt"]] <- PercentageFeatureSet(SR, pattern = "^MT-")
#head(MN@meta.data, 5) #查看Seurat object的方式
VlnPlot(SR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(SR, feature1 = "nFeature_RNA", feature2 = "percent.mt")

#plot1 <- FeatureScatter(MN, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(MN, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1, plot2))

SR_qc <- subset(SR, subset = nFeature_RNA < 12000 & nFeature_RNA > 2000 & percent.mt < 15)
# 141->124 这一步的threshold一开始设的太严格了，导致丢掉了很多的细胞，其实不需要这样严格的
#但是要把pool的几个样本去掉
SR_qc<-SR[,!colnames(SR) %in% c("TPM_SRR4047260", "TPM_SRR4047278","TPM_SRR4047295","TPM_SRR4047315","TPM_SRR4047390","TPM_SRR4047332","TPM_SRR4047331", "TPM_SRR4047323")] 
#14131 features across 120 samples within 1 assay 

############ 2. Normalizing the data (为找HVG作准备) #################
#“LogNormalize”:lognormalize”全局缩放的归一化方法，通过总表达值对每个细胞的基因表达值归一化，并将其乘以缩放因子(默认为10,000)，最后对结果进行对数变换
SR_normed<- NormalizeData(SR_qc, normalization.method = "LogNormalize", scale.factor = 10000)

############ 3. Highly variable gene ##################
SR_normed <- FindVariableFeatures(SR_normed, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SR_normed), 10)

# plot variable features with and without labels
VariableFeaturePlot(SR_normed)
LabelPoints(plot = plot1, points = top10,repel = TRUE,xnudge =0, ynudge=0)  #repel是为了字和字之间的间隔


################ 4. PCA ###################
#Scaling data before PCA, 使用线性变换(“scaling”)处理数据，过程：改变每个基因的表达，使细胞间的平均表达为0；缩放每个基因的表达，使细胞间的差异为1，这样高表达基因就不会影响后续分析
#所以在scale之后就不能再找highly variable gene了！！！
all.genes <- rownames(SR_normed)
SR_scaled <- ScaleData(SR_normed, features = all.genes,vars.to.regress = "percent.mt")  
#SR_scaled <- ScaleData(SR_normed, features = all.genes)
                       #如果regress out mt genes会发现基因表达量都不保序了
                       
#没有规定all.genes的话默认只对前2000个基因做scaling, vars.to.regress是要单独regress against each feature的变量,即在不考虑mt gene 带来的true biological variance时使用
SR_pca <- RunPCA(SR_scaled, features = VariableFeatures(object = SR_scaled))
#看一下组成每个PC的关键基因
library(ggplot2)
print(SR_pca[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(SR_pca, dims = 1:2, reduction = "pca")
#SR_pca@meta.data$celltype<-as.factor(SR_pca@meta.data$celltype)

#前面run过PCA的object才能umap
SR_umap <- RunUMAP(object = SR_pca, reduction = "pca", dims = 1:30)  #会慢一些
DimPlot(object = SR_umap, reduction = "umap",  
        repel = TRUE) #+ NoLegend()    #label = TRUE, group.by = "celltype",

DimHeatmap(SR_pca, dims = 1:3, cells = 50, balanced = TRUE)
# cells means plotting the extreme cells on both end of the spectrum.
# head(SR_umap@meta.data,4)

# How many components to include in clustering?
#通过permutate 确定真正重要的features
SR_JS <- JackStraw(SR_pca, num.replicate = 100)
SR_JS <- ScoreJackStraw(SR_JS, dims = 1:20)
JackStrawPlot(SR_JS, dims = 1:15)
ElbowPlot(SR_JS)  #或者直接找拐点

########### 5. Clustering ################
SR_pca <- FindNeighbors(SR_pca, dims = 1:10)
SR_pca <- FindClusters(SR_pca, resolution = 0.5)
head(Idents(SR_pca), 5)   #查看前5个细胞的cluster，的确分成了3群，但是有分错群的细胞

#找到分类错误的细胞，对比SR_pca和SR_scaled的分群情况
new.cluster.ids <- c("MN","iPSC","NPC")
names(new.cluster.ids) <- levels(SR_pca)
SR_pca <- RenameIdents(SR_pca, new.cluster.ids)
which(Idents(SR_pca)!=Idents(SR_scaled))

#subset 可以选择一部分cell和一部分feature
SR_pca<-subset(SR_pca, cells=which(Idents(SR_pca)==Idents(SR_scaled)))
#剩下116个细胞

#或者先用所有基因做一下也可以, 查看原始count矩阵
SR_pca@assays$RNA@counts  #自始至终不会变
SR_pca@assays$RNA@data   #在normalize的那一步改变，标准化数据的存放位置,没有改变基因表达的细胞间排序
SR_scaled@assays$RNA@scale.data  #scale之后出现，含负值

filtered<-gsub("TPM_", "", colnames(SR_pca@assays$RNA@scale.data))
write.table (filtered, file = "/Users/chenxinyi/Desktop/filtered116")
#保存seurat object为一个R variable
saveRDS(SR_pca, file = "/Users/chenxinyi/Desktop/filtered116.rds")

############### 6. correlation analysis ################

##太多了根本算不完
#把NA的细胞删除掉再做correlation
MN_expr<-SR_pca@assays$RNA@scale.data
MN_expr<-as.data.frame(t(MN_expr))  # 做转置时是matrix格式，需要再转化成dataframe一次
rownames(MN_expr)<-filtered
MN_se_psi= read.csv("/Users/chenxinyi/Desktop/MN_se_bi_psi.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
rownames(MN_se_psi)

ids<-rownames(MN_expr) 
MN_expr<-cbind(MN_expr,ID=ids)
ids<-rownames(MN_se_psi) 
MN_se_psi<-cbind(MN_se_psi,ID=ids)

merged<-merge(MN_expr,MN_se_psi,all=F,by="ID")
#48行，17124列=ID + 13674 gene + 3449 isoforms

cor_mat<-data.frame(matrix(NA,13674,3449))
rownames(cor_mat)<-colnames(merged[2:13675])
colnames(cor_mat)<-colnames(merged[13676:17124])

for (i in 2:13675){
  gene<-merged[,i]
  for (j in 13676:17124){
    as<-merged[,j]
    gene_na<-merged[which(!is.na(gene) & !is.na(as)),i]
    as_na<-merged[which(!is.na(gene) & !is.na(as)),j]
    cor_mat[i-1,j-13675]<-cor(gene_na,as_na,method = "spearman")
  }
}
