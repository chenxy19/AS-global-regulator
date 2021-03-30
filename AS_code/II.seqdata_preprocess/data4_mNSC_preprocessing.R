########### dataset 4 preprocessing ############

# reanalysing mouseSC2 with sce

mouseSC2_RSEM_nreads= read.csv("mouseSC2/mouseSC2_RSEM_nreads.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
colnames(mouseSC2_RSEM_nreads)<- gsub("expected_count_","",colnames(mouseSC2_RSEM_nreads))

ID = lapply(rownames(mouseSC2_RSEM_nreads),ID_convert)
rownames(mouseSC2_RSEM_nreads)<-ID
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

# merge three replicates for each single cell.
SRR_map=read.csv('mouseSC2/mouseSC2_descriptor_SRR.txt',sep="\t",stringsAsFactors = F)
SRR_map<-SRR_map[,c("sample_accession","run_accession","sample_title")]

title_convert<-function(ID){
  ID<-strsplit(ID, " ")[[1]][1]
  return(ID)
}

SRR_map$sample_title<-lapply(SRR_map$sample_title,title_convert)
SRR_map<-SRR_map[which(SRR_map$sample_title=="single_cell_RNAseq"),]
SAM_ID<-unique(SRR_map$sample_accession)

SAM_mat<-data.frame(matrix(NA,nrow(mouseSC2_RSEM_nreads),length(SAM_ID)))
colnames(SAM_mat)<-SAM_ID
rownames(SAM_mat)<-rownames(mouseSC2_RSEM_nreads)
for (samp in colnames(SAM_mat)){
  print(samp)
  SRR_id <- SRR_map[which(SRR_map$sample_accession==samp),"run_accession"]
  if (length(intersect(SRR_id,colnames(mouseSC2_RSEM_nreads)))!=0){
    SRR_id<-intersect(SRR_id,colnames(mouseSC2_RSEM_nreads))
    SRR_mat <- mouseSC2_RSEM_nreads[,SRR_id]
    if (length(SRR_mat)<10){
      SAM_mat[,samp]<-apply(SRR_mat,1,sum)
    }
    else{
      SAM_mat[,samp]<-SRR_mat
    }
  }
}

gene_ID<-unique(features[rownames(SAM_mat),"external_gene_name"])

gene_ID <- gene_ID[!is.na(gene_ID)]
SAM_mat2<-data.frame(matrix(NA,length(gene_ID),ncol(SAM_mat)))
colnames(SAM_mat2)<-colnames(SAM_mat)
rownames(SAM_mat2)<-gene_ID
for (samp in gene_ID){
  SRR_id <- features[which(features$external_gene_name==samp),"ensembl_gene_id"]
  if (length(SRR_id)==1){
    SAM_mat2[samp,]<-SAM_mat[SRR_id,]
  }else{
    if (length(intersect(SRR_id,rownames(SAM_mat)))!=0){
      SRR_id<-intersect(SRR_id,rownames(SAM_mat))
      SRR_mat <- SAM_mat[SRR_id,]
      SAM_mat2[samp,]<-apply(SRR_mat,2,sum)
      print(samp)
    }
  }
}

SAM_mat2<-rbind(SAM_mat2, SAM_mat[grepl("ERCC",rownames(SAM_mat)),])

write.table(SAM_mat2, file="mouseSC2/mouseSC2_nreads_processed.csv")

# no pheno annotation









is.spike <- grepl("^ERCC-", gene.names)
summary(is.spike)
is.mito <- !is.na(features[!is.spike,]$chromosome_name) & features[!is.spike,]$chromosome_name=="MT"
summary(is.mito)


library(SingleCellExperiment)
sce <- SingleCellExperiment(list(counts=as.matrix(SAM_mat)))
sce <- splitAltExps(sce, ifelse(is.spike, "ERCC", "gene"))
sce <- splitAltExps(sce, ifelse(is.mito, "mito", "gene"))
sce
#55364 150 old cells
colData(sce) <- perCellQCMetrics(sce) 
colnames(colData(sce))

rowData(sce)$entrezID<-rownames(sce)
rowData(sce)$symbol<-lapply(rownames(sce),map_name)
rownames(sce)<-rowData(sce)$symbol
colData(sce)$Gapdh<-counts(sce)["Gapdh",]
colData(sce)$Actb<-counts(sce)["Actb",]
colData(sce)$Slc1a3<-counts(sce)["Slc1a3",]

multiplot(cols=3,
          plotColData(sce,y="sum"),   # sum of counts, or library size
          plotColData(sce,y="detected"),  # number of detected features
          plotColData(sce,y="altexps_mito_percent"),
          plotColData(sce,y="Gapdh"),
          plotColData(sce,y="Actb"),
          plotColData(sce,y="Slc1a3")
)

# also remove cells with high numbers of genes detected (possibly doublets)
libsize.drop <- isOutlier(sce$sum, nmads=2, type="lower", log=TRUE) 
feature.drop <- isOutlier(sce$detected, nmads=2, type="lower", log=TRUE) 
#mito.drop <- isOutlier(sce$altexps_mito_percent, nmads=3, type="higher")
gapdh.drop <- isOutlier(sce$Gapdh, nmads=2, type="lower")
actb.drop <- isOutlier(sce$Actb, nmads=2, type="lower")
Slc1a3.drop <- isOutlier(sce$Slc1a3, nmads=3, type="higher")

discard <- libsize.drop | feature.drop | gapdh.drop | actb.drop | Slc1a3.drop

data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), Bygapdh=sum(gapdh.drop),
           Byactb=sum(actb.drop),BySlc1a3=sum(Slc1a3.drop), Remaining.in.batch=sum(!discard))
# ByLibSize ByFeature Bygapdh Byactb BySlc1a3 Remaining.in.batch
# 1        16         5       0      0        2                128

sce <- sce[, !discard]
dim(sce) 
# 55364   128

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
# keep 32675 genes

ncells <- scater::nexprs(sce, byrow=TRUE)
rowData(sce)$Ncells <- ncells
plot(ave.counts[keep], ncells, xlab="Average count", ylab="Number of cells", log="x", pch=16, cex=0.5)

#the size factors are applied to compute normalized log-expression values.
set.seed(10000)
clusters <- quickCluster(sce, method="igraph", min.mean=1)
sce <- computeSumFactors(sce, cluster=clusters, min.mean=1)
summary(sizeFactors(sce))

# correlates with library size
plot(sizeFactors(sce), sce$sum/1e6, log="xy", ylab="Library size (millions)", 
     xlab="Size factor")
sce <- logNormCounts(sce,log = TRUE)
saveRDS(sce, file="mouseSC2/sce_all.rds")

# PCA
sce<- denoisePCA(sce, var.out, assay.type="logcounts")
ncol(reducedDim(sce))
pc.out <- attr(reducedDim(sce), "percentVar")
plot(pc.out)

set.seed(100)
colData(sce)$qNSC1m<-logcounts(sce["Cxcl14",])
colData(sce)$qNSC2m<-logcounts(sce["Cpe",])
colData(sce)$aNSC1m<-logcounts(sce["Rplp1",])
colData(sce)$aNSC2m<-logcounts(sce["Dek",])

sce <- runTSNE(sce, use_dimred="PCA", perplexity=30)
pca1 <- plotPCA(sce, colour_by="qNSC1m") + fontsize  
pca2 <- plotPCA(sce, colour_by="qNSC2m") + fontsize  
pca3 <- plotPCA(sce, colour_by="aNSC1m") + fontsize  
pca4 <- plotPCA(sce, colour_by="aNSC2m") + fontsize  
multiplot( pca1,pca2,pca3,pca4, cols=2)

# SC3 clustering
rowData(sce)$feature_symbol<-rowData(sce)$symbol
sce<-sc3(sce,ks=4)

# plot known differential genes in 4 clusters to confirm
gene2conv<-c("ENSMUSG00000021508","ENSMUSG00000036949","ENSMUSG00000024810","ENSMUSG00000026385","ENSMUSG00000024411","ENSMUSG00000067274","ENSMUSG00000037894","ENSMUSG00000028832")
shuffled.sce <- sce[, sample(ncol(sce))]
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
scater::plotExpression(shuffled.sce, features[gene2conv,"external_gene_name"], 
                       colour_by="sc3_4_clusters", show_violin = FALSE) + fontsize

part<-rownames(colData(sce))[which(colData(sce)$sc3_4_clusters==3)]
sce_part<-sce[,part]

#convert to human name, 去掉重复基因
rowData(sce_part)<-cbind(rowData(sce_part), features[rowData(sce_part)$entrezID,])
sce_part<-sce_part[which(!is.na(rowData(sce_part)$HGNC.symbol) & !duplicated(rowData(sce_part)$HGNC.symbol)),]
rownames(sce_part)<-rowData(sce_part)$HGNC.symbol
#12938 76

norm_expr<-logcounts(sce_part)
write.table(norm_expr, file="mouseSC2/mouseSC2_sce_expr_qNSC1.csv")
