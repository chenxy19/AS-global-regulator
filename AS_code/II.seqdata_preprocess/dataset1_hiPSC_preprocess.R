# sce dataset1 preprocessing

expr_mat= read.csv("iPSC/iPSC_gene_nreads.csv",header = TRUE,stringsAsFactors = FALSE) 
colnames(expr_mat)<- gsub("NumReads_","",colnames(expr_mat))
colnames(expr_mat)[2]<-"SRR4047457" 

#add up transcripts into genes (每一列都求和)
gene_nreads <- group_by(expr_mat, X) %>% summarize_all(sum)
gene_nreads<-as.data.frame(gene_nreads)
rownames(gene_nreads)<-gene_nreads$X
gene_nreads<-gene_nreads[,-1]
write.table(gene_nreads,file="iPSC/iPSC_gene_nreads.csv")

# merge all 3 cell types

nreads1= read.csv("iPSC/MN_nreads_processed.csv",sep=" ",header = TRUE,stringsAsFactors = FALSE) 
nreads2= read.csv("iPSC/NPC_nreads_processed.csv",sep=" ",header = TRUE,stringsAsFactors = FALSE) 
nreads3= read.csv("iPSC/iPSC_nreads_processed.csv",sep=" ",header = TRUE,stringsAsFactors = FALSE) 
nreads4<-merge_by_rownames(nreads1,nreads2)
nreads5<-merge_by_rownames(nreads4,nreads3)
write.table(nreads5,file="iPSC/all_nreads_processed.csv")
colnames(nreads)
 
nreads6<-as.data.frame(t(nreads5))
nreads6$cell_type<-0
nreads6$cell_type[1:ncol(nreads1)]<-"MN"
nreads6$cell_type[(ncol(nreads1)+1):ncol(nreads4)]<-"NPC"
nreads6$cell_type[(ncol(nreads4)+1):ncol(nreads5)]<-"iPSC"
nreads6<-data.frame(nreads6$cell_type,row.names = rownames(nreads6))
colnames(nreads6)<-"cell_type"
write.table(nreads6,file="iPSC/all_pheno.csv")
