########## data3_mHSC preprocessing ########

"""
SCE analysis for the mouseSC data (data3) without ERCC spike-in
repeat the SC3 clustering in the article, and identify 4month WT cells as a pure population.
"""

# !!! Note !!!
# An easy method to account for missing and duplicate symbols. The code below will replace missing symbols with the Ensembl identifier and concatenate duplicated symbols with the (unique) Ensembl identifiers.
#library(scater)
#rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ENSEMBL, rowData(sce)$SYMBOL)
#head(rownames(sce))

mouseSC_RSEM_nreads= read.csv("mouseSC/mouseSC_RSEM_nreads.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
colnames(mouseSC_RSEM_nreads)<- gsub("expected_count_","S",colnames(mouseSC_RSEM_nreads))

ID_convert<-function(ID){
  ID<-strsplit(ID, "[.]")[[1]][1]
  return(ID)
}

ID = lapply(rownames(mouseSC_RSEM_nreads),ID_convert)
rownames(mouseSC_RSEM_nreads)<-ID
gene.names<-as.vector(unlist(ID))

library(biomaRt)
# To know which BioMart databases are available see the listMarts function. To know which datasets are available within a BioMart database, first select the BioMart database using useMart and then use the listDatasets function on the selected BioMart, see listDatasets function.
ensembl <- useMart(biomart='ensembl', dataset="mmusculus_gene_ensembl")

ensemblGenes <- getBM(attributes=c('ensembl_gene_id',  'chromosome_name', 
                                   'external_gene_name'), filters="ensembl_gene_id", 
                      values=gene.names, mart=ensembl) 

human = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl",dataset = "mmusculus_gene_ensembl")
#genes=c("Zfp286","Tmx2")
homo_genes=getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
                  values = ensemblGenes$external_gene_name, mart=mouse,attributesL = c("hgnc_symbol","ensembl_gene_id"),
                  martL = human,
                  uniqueRows = T)
homo_genes_dedup<- homo_genes %>% distinct(MGI.symbol, .keep_all = TRUE)
                                           
features <- ensemblGenes[match(gene.names, ensemblGenes$ensembl_gene_id),]
features$ensembl_gene_id <- gene.names
row.names(features) <- gene.names
features<-merge(features,homo_genes_dedup, by.x="external_gene_name",by.y="MGI.symbol", all=TRUE)
rownames(features)<-features$ensembl_gene_id
#把小鼠基因名和同源基因的基因名都放进去,如果小鼠基因map到多个人基因上则只保留一个

gene_ID<-unique(features[rownames(mouseSC_RSEM_nreads),"external_gene_name"])

gene_ID <- gene_ID[!is.na(gene_ID)]
SAM_mat<-data.frame(matrix(NA,length(gene_ID),ncol(mouseSC_RSEM_nreads)))
colnames(SAM_mat)<-colnames(mouseSC_RSEM_nreads)
rownames(SAM_mat)<-gene_ID
for (samp in gene_ID){
  SRR_id <- features[which(features$external_gene_name==samp),"ensembl_gene_id"]
  if (length(SRR_id)==1){
    SAM_mat[samp,]<-mouseSC_RSEM_nreads[SRR_id,]
  }else{
    if (length(intersect(SRR_id,rownames(mouseSC_RSEM_nreads)))!=0){
      SRR_id<-intersect(SRR_id,rownames(mouseSC_RSEM_nreads))
      SRR_mat <- mouseSC_RSEM_nreads[SRR_id,]
      SAM_mat[samp,]<-apply(SRR_mat,2,sum)
      print(samp)
    }
  }
}

write.table(SAM_mat, file="mouseSC/mouseSC_nreads_processed.csv")

SRR_map=read.csv('mouseSC/mouseSC_descriptor_extra.txt',sep=",",stringsAsFactors = F)
SRR_map<-SRR_map[,c("Run","AGE","Genotype")]
rownames(SRR_map)<-SRR_map$Run

SRR_map$cell_type<-0
SRR_map$cell_type[which(SRR_map$AGE=="4 Months (young)" & SRR_map$Genotype=="Wild Type")]<-"4W"
SRR_map$cell_type[which(SRR_map$AGE=="12 Months" & SRR_map$Genotype=="Wild Type")]<-"12W"
SRR_map$cell_type[which(SRR_map$AGE=="18 Months (Old)" & SRR_map$Genotype=="Wild Type")]<-"18W"
SRR_map$cell_type[which(SRR_map$AGE=="4 Months (young)" & SRR_map$Genotype=="JAK2 V617F Homozygous")]<-"4M"
SRR_map$cell_type[which(SRR_map$AGE=="12 Months"  & SRR_map$Genotype=="JAK2 V617F Homozygous")]<-"12M"
SRR_map$cell_type[which(SRR_map$AGE=="18 Months (Old)" & SRR_map$Genotype=="JAK2 V617F Homozygous")]<-"18M"

write.table(SRR_map, file="mouseSC/mouseSC_pheno.csv")





#clus1<-read.csv("mouseSC/mouseSC_clus1_expr.csv",sep = " ",stringsAsFactors = F)
#rownames(SRR_map)<-SRR_map$Run
#length(SRR_map[colnames(sce)[colData(sce)$seurat_cluster=="0"],"BioSample"])
# Seurat cluster 0 is basically the clus1 we analyzed before



is.spike <- grepl("^ERCC-", gene.names)
summary(is.spike)
#92个spike in (确实有ERCC开头的真基因)
is.mito <- !is.na(features[!is.spike,]$chromosome_name) & features[!is.spike,]$chromosome_name=="MT"
summary(is.mito)
#37个mitochondria gene

library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=as.matrix(mouseSC_RSEM_nreads)))
sce <- splitAltExps(sce, ifelse(is.spike, "ERCC", "gene"))
sce <- splitAltExps(sce, ifelse(is.mito, "mito", "gene"))
sce
#55364 1139
colData(sce) <- perCellQCMetrics(sce) 
colnames(colData(sce))


multiplot(cols=2,
          plotColData(sce,  x="group",y="sum"),   # sum of counts, or library size
          plotColData(sce,  x="group",y="detected"),  # number of detected features
          plotColData(sce,  x="group",y="altexps_ERCC_percent"),
          plotColData(sce,  x="group",y="altexps_mito_percent")
)

#原文的QC太粗暴了，用outlier的方法来排除

libsize.drop <- isOutlier(sce$sum, nmads=3, type="lower", log=TRUE, batch=sce$group) 
feature.drop <- isOutlier(sce$detected, nmads=3, type="lower", log=TRUE, batch=sce$group) 
mito.drop <- isOutlier(sce$altexps_mito_percent, nmads=3, type="higher", batch=sce$group)
spike.drop <- isOutlier(sce$altexps_ERCC_percent, nmads=3, type="higher", batch=sce$group)
discard <- libsize.drop | feature.drop | spike.drop | mito.drop
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), ByMito=sum(mito.drop),
           BySpike=sum(spike.drop), Remaining.in.batch=sum(!discard))
sce <- sce[, !discard]
dim(sce) 

rowData(sce)$entrezID<-rownames(sce)
map_name<-function(id){
  symbol<-features$external_gene_name[which(features$ensembl_gene_id==id)] 
  return(symbol)
}
rowData(sce)$symbol<-lapply(rownames(sce),map_name)
#比写循环快得多
rownames(sce)<-rowData(sce)$symbol

# check if highly expressed genes are correct
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
scater::plotHighestExprs(sce, n=50) + fontsize

# plot mean expression level
ave.counts <- apply(counts(sce),1,mean)
rowData(sce)$AveCount <- ave.counts
hist(log10(ave.counts), breaks=100, main="", col="grey", xlab=expression(Log[10]~"average count"))

# filter gene not expressed
keep <- ave.counts > 0
sce <- sce[keep,]
summary(keep)

ncells <- scater::nexprs(sce, byrow=TRUE)
rowData(sce)$Ncells <- ncells
plot(ave.counts, ncells, xlab="Average count", ylab="Number of cells", log="x", pch=16, cex=0.5)

#the size factors are applied to compute normalized log-expression values.
set.seed(10000)
clusters <- quickCluster(sce, method="igraph", min.mean=1)
sce <- computeSumFactors(sce, cluster=clusters, min.mean=1)
summary(sizeFactors(sce))

# correlates with library size
plot(sizeFactors(sce), sce$sum/1e6, log="xy", ylab="Library size (millions)", 
     xlab="Size factor", col=ifelse(sce$phenotype=="naive", "black", "grey"))

sce_ERCC <- computeSumFactors(altExp(sce,"ERCC"), cluster=clusters, min.mean=1)
#computeSpikeFactors(sce, type="ERCC", general.use=FALSE) 
plot(sizeFactors(sce), sizeFactors(sce_ERCC), log="xy", ylab="Spike-in size factor", 
     xlab="Deconvolution size factor")

sce <- logNormCounts(sce,log = TRUE)
saveRDS(sce, file="mouseSC/sce_all.rds")
#sce <- readRDS("results-preprocess/sce_all.rds")

# for PCA and clustering
# 12月和4月+18月之间的batch effect特明显,但是因为technical variation和biological variation混在一起，这样其实特不好
#scater::norm_exprs(sce) <- removeBatchEffect(logcounts(sce),
#                design=model.matrix(~sce$group), batch=sce$AGE)
#assayNames(sce) 

#单独给同一批测序的4W和18W聚类
part<-rownames(colData(sce))[which(colData(sce)$group=="4W" | colData(sce)$group=="18W")]
sce_part<-sce[,part]
# 24005 280
sce_part<- denoisePCA(sce_part, var.out, assay.type="logcounts")
ncol(reducedDim(sce_part))
pc.out <- attr(reducedDim(sce_part), "percentVar")
plot(pc.out)

set.seed(100)
sce_part <- runTSNE(sce_part, use_dimred="PCA", perplexity=30)
tsne1 <- scater::plotTSNE(sce_part, colour_by="AGE") +  fontsize
pca1 <- plotPCA(sce_part, colour_by="AGE") + fontsize
#tsne3 <- scater::plotTSNE(sce, colour_by="batch") + fontsize
#tsne4 <- scater::plotTSNE(sce, colour_by="phenotype") + fontsize
multiplot(tsne1, pca1, cols=2)

# SC3 clustering
rowData(sce_part)$feature_symbol<-rowData(sce_part)$symbol
sce_part<-sc3(sce_part,ks=3:5)

sum(sce_part$sc3_5_clusters==2 & sce_part$group=="4W")-sum(sce_part$sc3_5_clusters==2)
#第二个cluster 确实是18W specific cluster

#convert to human name, 去掉重复基因
rowData(sce_part)<-cbind(rowData(sce_part), features[rowData(sce_part)$entrezID,])
sce_part<-sce_part[which(!is.na(rowData(sce_part)$HGNC.symbol) & !duplicated(rowData(sce_part)$HGNC.symbol)),]
rownames(sce_part)<-rowData(sce_part)$HGNC.symbol
#13418 280 

#12month的数据ERCC含量比较高，而且和4month, 18month明显有区分，就不用这部分data了。
#所以说按cluster分的意义其实不大，因为4month和18month都可以分成更细的小群。还是对4 month和18month的细胞做一下rank
sce_4W <- sce_part[,which(colData(sce_part)$group=="4W")]
sce_18W <- sce_part[,which(colData(sce_part)$group=="18W")]
norm_expr_4W<-logcounts(sce_4W)
norm_expr_18W<-logcounts(sce_18W)
write.table(norm_expr_4W, file="mouseSC/mouseSC_norm_expr_4W.csv")
write.table(norm_expr_18W, file="mouseSC/mouseSC_norm_expr_18W.csv")



