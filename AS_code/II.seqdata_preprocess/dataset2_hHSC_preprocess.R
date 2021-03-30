########## data2_hHSC preprocessing ########

#换用kallisto求表达矩阵

nreads= read.csv("HSC/hHSC_kal_nreads.csv",header = TRUE,stringsAsFactors = FALSE) 
colnames(nreads) <- gsub("est_counts_","",colnames(nreads))
gene_nreads <- group_by(nreads, X) %>% summarize_all(sum) # takes a long time
gene_nreads<-as.data.frame(gene_nreads)
rownames(gene_nreads)<-gene_nreads$X
gene_nreads<-gene_nreads[,-1]

# merge double runs into single cell data
SRR_map=read.csv('HSC/HSC_descriptor_SRR.txt',sep="\t",stringsAsFactors = F)
SRR_map<-SRR_map[,c("sample_accession","run_accession","sample_title")]

SAM_ID<-unique(SRR_map$sample_accession)
SAM_mat<-data.frame(matrix(NA,nrow(gene_nreads),length(SAM_ID)))
colnames(SAM_mat)<-SAM_ID
rownames(SAM_mat)<-rownames(gene_nreads)
for (samp in colnames(SAM_mat)){
  SRR_id <- SRR_map[which(SRR_map$sample_accession==samp),"run_accession"]
  if (length(intersect(SRR_id,colnames(gene_nreads)))!=0){
    SRR_id<-intersect(SRR_id,colnames(gene_nreads))
    SRR_mat <- gene_nreads[,SRR_id]
    if (length(SRR_mat)<10){
      SAM_mat[,samp]<-apply(SRR_mat,1,sum)
      print(samp)
    }
    else{
      SAM_mat[,samp]<-SRR_mat
    }
  }
}
SAM_mat_nonna<-SAM_mat[, apply(SAM_mat, 2, function(y) any(!is.na(y)))]
# 192个细胞filter过剩下185个

# remove the version number in EnID
nreads<-SAM_mat_nonna
for (i in 1:nrow(nreads)){
  nreads$gene_name[i]<-strsplit(rownames(nreads)[i],split="[.]")[[1]][1]
}
# 有个别重复的基因删除掉
nreads_dedup<- nreads %>% distinct(gene_name, .keep_all = TRUE)
rownames(nreads_dedup) <- nreads_dedup$gene_name
nreads_dedup<-nreads_dedup[,-ncol(nreads_dedup)]
gene.names<-rownames(nreads_dedup)

# convert gene ID to name
library(BiocFileCache)
bmcache <- BiocFileCache("biomart", ask = FALSE)
loc <- bfcquery(bmcache, "hg38.ensGene", exact=TRUE)$rpath
if (length(loc)==0L) {
  library(biomaRt)
  ensembl <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset="hsapiens_gene_ensembl", 
                     host="aug2017.archive.ensembl.org") # Ensembl version 90.
  ensemblGenes <- getBM(attributes=c('ensembl_gene_id',  'chromosome_name', 'gene_biotype', 
                                     'external_gene_name', 'entrezgene'), filters="ensembl_gene_id", 
                        values=gene.names, mart=ensembl) 
  saveRDS(ensemblGenes, file=bfcnew(bmcache, "hg38.ensGene"))
} else {
  ensemblGenes <- readRDS(loc)
}

features <- ensemblGenes[match(gene.names, ensemblGenes$ensembl_gene_id),]
features$ensembl_gene_id <- gene.names
features$entrezgene <- gsub(" ", "", as.character(features$entrezgene))
row.names(features) <- gene.names
head(features)

#又有好多map到同一个基因上的,需要求和
gene_name_dedup<-features[rownames(nreads_dedup),"external_gene_name"]

gene_ID <- unique(gene_name_dedup)
gene_ID <- gene_ID[!is.na(gene_ID)]
SAM_mat<-data.frame(matrix(NA,length(gene_ID),ncol(nreads_dedup)))
colnames(SAM_mat)<-colnames(nreads_dedup)
rownames(SAM_mat)<-gene_ID
for (samp in gene_ID){
  SRR_id <- features[which(features$external_gene_name==samp),"ensembl_gene_id"]
  if (length(SRR_id)==1){
    SAM_mat[samp,]<-nreads_dedup[SRR_id,]
  }else{
    if (length(intersect(SRR_id,rownames(nreads_dedup)))!=0){
      SRR_id<-intersect(SRR_id,rownames(nreads_dedup))
      SRR_mat <- nreads_dedup[SRR_id,]
      SAM_mat[samp,]<-apply(SRR_mat,2,sum)
      print(samp)
  }
  }
}
# prevent the loss of ERCC genes
SAM_mat<-rbind(SAM_mat,nreads_dedup[grepl("^ERCC",rownames(nreads_dedup)),])
write.table(SAM_mat, file="HSC/hHSC_nreads_processed.csv")


SRR_map$cell_type<-as.factor(1+as.numeric(grepl("patient 2",SRR_map$sample_title)))
SRR_map<-SRR_map[,c("sample_accession","cell_type")]
SRR_map<-unique(SRR_map)
rownames(SRR_map)<-SRR_map$sample_accession
SRR_map<-SRR_map[which(rownames(SRR_map) %in% colnames(SAM_mat_nonna)),]
write.table(SRR_map, file="HSC/hHSC_pheno.csv")









# identify spike gene and mito gene
is.spike <- grepl("^ERCC-", gene.names)
summary(is.spike)
#92个spike in (确实有ERCC开头的真基因)

is.mito <- !is.na(features[!is.spike,]$chromosome_name) & features[!is.spike,]$chromosome_name=="MT"
summary(is.mito)
#37个mitochondria gene,折腾了半天，其实和grep "^MT-"找出来的一样
# 为啥kallisto弄完了连个mitochondria 基因都没有？？？

library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=as.matrix(nreads)))
sce <- splitAltExps(sce, ifelse(is.spike, "ERCC", "gene"))
sce <- splitAltExps(sce, ifelse(is.mito, "mito", "gene"))

sce
#60660 185

#QC
colData(sce) <- perCellQCMetrics(sce) 
colnames(colData(sce))

p_map<-unique(SRR_map[,c("sample_accession","sample_title")])
p1<-p_map$sample_accession[grepl("patient 1",p_map$sample_title)]
p2<-p_map$sample_accession[grepl("patient 2",p_map$sample_title)]
colData(sce)$patient <- as.factor(1 + as.numeric(rownames(colData(sce)) %in% p2))

multiplot(cols=2,
          plotColData(sce,  x="patient",y="sum"),   # sum of counts, or library size
          plotColData(sce,  x="patient",y="detected"),  # number of detected features
          plotColData(sce,  x="patient",y="altexps_ERCC_percent")
          #plotColData(sce,  x="patient",y="altexps_mito_percent")
)

feature.drop.p1<- colData(sce)$patient == "1" & colData(sce)$detected <2000
#45 个
feature.drop.p2<- colData(sce)$patient == "2" & colData(sce)$detected <3500
#8个
discard<-feature.drop.p1 | feature.drop.p2
sce <- sce[, !discard]
ERCC.mito.drop<-colData(sce)$altexps_ERCC_percent > 20 
sce <- sce[, !ERCC.mito.drop]
#132个cell, 44个patient1, 88个patient 2


rowData(sce)$entrezID<-rownames(sce)
map_name<-function(id){
  symbol<-features$external_gene_name[which(features$ensembl_gene_id==id)] 
  return(symbol)
}
rowData(sce)$symbol<-lapply(rownames(sce),map_name)
#比写循环快得多
rownames(sce)<-rowData(sce)$symbol

fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
scater::plotHighestExprs(sce, n=50) + fontsize

sce<-sce[!is.na(rownames(sce)),]
#  57158 132 
ave.counts <- apply(counts(sce),1,mean)
rowData(sce)$AveCount <- ave.counts
hist(log10(ave.counts), breaks=100, main="", col="grey", xlab=expression(Log[10]~"average count"))

keep <- ave.counts > 0
sce <- sce[keep,]
summary(keep)
# 31788 left

#size factor normalization to account for seq depth variability
set.seed(10000)
clusters <- quickCluster(sce, method="igraph", min.mean=1)
sce <- computeSumFactors(sce, cluster=clusters, min.mean=1)
summary(sizeFactors(sce))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.05752 0.67198 1.01563 1.00000 1.26477 2.18881 

sce2 <- computeSpikeFactors(sce, "ERCC")
summary(sizeFactors(sce2))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.01405  0.35781  0.53762  1.00000  0.83027 18.49483

sce <- logNormCounts(sce,size_factors=sizeFactors(sce)*sizeFactors(sce2),log =TRUE)

sce<- denoisePCA(sce, var.out, assay.type="logcounts")
ncol(reducedDim(sce))
pc.out <- attr(reducedDim(sce), "percentVar")
plot(pc.out)

set.seed(100)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=30)
tsne1 <- scater::plotTSNE(sce, colour_by="patient") +  fontsize
pca1 <- plotPCA(sce, colour_by="patient") + fontsize
tsne2 <- scater::plotTSNE(sce, colour_by="sc3_3_clusters") + fontsize
#tsne4 <- scater::plotTSNE(sce, colour_by="phenotype") + fontsize
multiplot(tsne1, tsne2, cols=2)

rowData(sce)$feature_symbol<-rowData(sce)$symbol
sce<-sc3(sce,ks=2:3)

colLabels(sce) <- colData(sce)$sc3_3_clusters
markers <- findMarkers(sce)
# marker gene 跟原文完全不一样。。。
