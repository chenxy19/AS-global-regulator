"""
1. Aggregate iPSC expression matrix and mean psi
   input: expr_mat (from single_cell_pathway_activity.R); all_mean_psi_corrected0831.csv
   output: iPSC_expr_psi.csv

2. Cell cycle analysis (Seurat assign cell cycle to cells based on marker genes)
   input: Seurat object recontructed from gene_mat2.csv
   
3. Convert human cell cycle genes to mouse cell cycle genes 
  
"""

############ 1. Aggregate iPSC expression matrix and mean psi ##########
expr_mat<-t(expr_mat)
rownames(expr_mat)<-gsub("TPM_","",rownames(expr_mat))
all_psi_new= read.csv("/Users/chenxinyi/Desktop/output/all_mean_psi_corrected0831.csv",sep=" ",header = TRUE,stringsAsFactors = FALSE) 
#坏处是这里只有111个基因，有时间可以重新把别的细胞也求一下
all_expr_psi<-merge_by_rownames(expr_mat,all_psi_new)
iPSC_expr_psi<-all_expr_psi[which(all_expr_psi$cell_type=="iPSC"),]
write.table (t(iPSC_expr_psi), file = "/Users/chenxinyi/Desktop/output/iPSC_expr_psi.csv")


############ 2. cell cycle analysis ##########
iPSC_expr_psi= read.csv("/Users/chenxinyi/Desktop/output/iPSC_expr_psi.csv",sep=" ",header = TRUE,stringsAsFactors = FALSE) 
# 13677 * 33

gene_mat2= read.csv("/Users/chenxinyi/Desktop/output/gene_mat2.csv",sep=" ",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 

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


############ 2.1 QC of cells ################
#the set of all genes starting with MT- as a set of mitochondrial genes
SR[["percent.mt"]] <- PercentageFeatureSet(SR, pattern = "^MT-")
SR[["sample_name"]] <- gsub("^TPM_","",colnames(SR))
#head(MN@meta.data, 5) #查看Seurat object的方式
VlnPlot(SR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(SR, feature1 = "nFeature_RNA", feature2 = "percent.mt")

SR_qc <- subset(SR, subset = nFeature_RNA < 12000 & nFeature_RNA > 2000 & percent.mt < 15)
SR_qc2 <- subset(SR_qc, subset= sample_name %in% colnames(iPSC_expr_psi))

############ 2.2 Normalizing the data (为找HVG作准备) #################
#“LogNormalize”:lognormalize”全局缩放的归一化方法，通过总表达值对每个细胞的基因表达值归一化，并将其乘以缩放因子(默认为10,000)，最后对结果进行对数变换
SR_normed<- NormalizeData(SR_qc, normalization.method = "LogNormalize", scale.factor = 10000)

############ 2.3 Highly variable gene and cell cycle assignment ##################
SR_normed <- FindVariableFeatures(SR_normed, selection.method = "vst", nfeatures = 2000)
SR_normed <- ScaleData(SR_normed, features = rownames(SR_normed))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
SR_normed <- CellCycleScoring(SR_normed , s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(SR_normed[[]])
# Visualize the distribution of cell cycle markers across
RidgePlot(SR_normed, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase
PCA_genes<-setdiff(c(s.genes, g2m.genes),c("UHRF1", "MLF1IP", "CASP8AP2", "FAM64A", "HN1", "CENPA"))
SR_normed <- RunPCA(SR_normed, features = PCA_genes)
#很蠢的是只有细胞数多于feature数才能做PCA。所以只能把所有细胞都分析了。正好看看是不是成熟细胞在G1的少一些
DimPlot(SR_normed,reduction="pca")

#画iPSC 各个细胞周期的 psi的histogram
#sum(SR_normed[[]]$old.ident=="iPSC" & SR_normed[[]]$Phase=="S")
iPSC_normed <- subset(SR_normed, subset= old.ident=="iPSC")
iPSC_cc<-iPSC_normed[[]]
rownames(iPSC_cc)<-iPSC_cc$sample_name
iPSC_expr_psi<-as.data.frame(t(iPSC_expr_psi))
iPSC_merged <- merge_by_rownames(iPSC_cc, iPSC_expr_psi)

hist(as.numeric(iPSC_merged$glbl_mean_psi[which(iPSC_merged$Phase=="S")]),
     col=rgb(1,0,0,1/4),xlab = "global mean psi", 
     main="psi in different cell cycle stages in iPSC",breaks = 20)
hist(as.numeric(iPSC_merged$glbl_mean_psi[which(iPSC_merged$Phase=="G2M")]),
     add=T,col=rgb(0,0,1,1/4),xlab = "global mean psi")
legend("topright", c("S", "G2M"),fill=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)),bty="n")  
#bty="n"是legend周围的边框不要

# all three types of cells
all_SR_normed<-SR_normed[[]]
rownames(all_SR_normed) <- gsub("^TPM_","",rownames(all_SR_normed))
all_phase_psi<-merge_by_rownames(all_SR_normed,all_psi_new)

cell_type="iPSC"
hist(as.numeric(all_phase_psi$glbl_mean_psi[which(all_phase_psi$Phase=="G1")]),
     col=rgb(1,0,0,1/4),xlab = "global mean psi", xlim = range(0.65,0.75),
     main="psi in different cell cycle stages in all cells",breaks = 20)
hist(as.numeric(all_phase_psi$glbl_mean_psi[which(all_phase_psi$Phase=="G2M")]),
     add=T,main=paste("psi in different cell cycle stages in",cell_type,sep=" "),
     col=rgb(0,0,1,1/4),xlab = "global mean psi",breaks = 20)
hist(as.numeric(all_phase_psi$glbl_mean_psi[which(all_phase_psi$Phase=="S")]),
     add=T,col=rgb(0,1,0,1/4),xlab = "global mean psi",breaks = 20)
legend("topleft", c("G1","G2M","S"),fill=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4),rgb(0,1,0,1/4)),bty="n")  

G1<-as.numeric(all_phase_psi$glbl_mean_psi[which(all_phase_psi$Phase=="G1" & all_phase_psi$cell_type=="NPC")]) 
G2M<-as.numeric(all_phase_psi$glbl_mean_psi[which(all_phase_psi$Phase=="G2M" & all_phase_psi$cell_type=="NPC")]) 
S<-as.numeric(all_phase_psi$glbl_mean_psi[which(all_phase_psi$Phase=="S" & all_phase_psi$cell_type=="NPC")]) 

shapiro.test(G1)
shapiro.test(G2M)
shapiro.test(S)
#大于0.05 则无法拒绝其为正态分布
#bartlett test是方差齐性检验
t.test(G1,G2M,paired = FALSE,var.equal = F)  #没有显著区别


########### 3. human cell cycle marker genes converted to mice #############
library(biomaRt)
#目前ensembl网站在修缮。。。should be able to use before...
human = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
class(human)
mouse = useMart("ensembl",dataset = "mmusculus_gene_ensembl")
class(mouse)

cc_genes <- cc.genes$s.genes
#listAttributes函数可以检索可能的属性列表，listFilters可以检索可能的过滤器列表
cc_mgene <- getLDS(attributes = c("hgnc_symbol"), filters="hgnc_symbol",
                   values=genes, mart=human,
                   attributesL = c("mgi_symbol","chromosome_name","start_position"),
                   martL=mouse,
                   uniqueRows=T)

#只好勉为其难的转化成小写试试
cc_mgenes<-lapply(cc_genes,tolower)
