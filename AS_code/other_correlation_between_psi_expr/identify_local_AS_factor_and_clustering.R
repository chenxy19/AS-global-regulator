"""
find local regulation factors first by clustering AS events.

1. clustering to find local factors: only take AS that are detected in >=80 cells, and fill the NAs by sampling from a background distribution.
   Then calculate pairwise distance between AS (Euclidean distance).
   Input: MN_psi/se/psi.csv, the psi for each AS in each cell.
   Output: all_se_psi_filt100cells.csv, AS_distance.csv, all_se_psi_nonna.csv.
   
1.5 clustering only AS in iPSC

2. Find GO for different clusters of AS
   Input: AS80_label.csv, from clustering result on Python.
   
3. Use mean psi for a cluster of AS to identify local factors.
   Input: event_label.csv
   Output: grp1_list, grp2_list, NPC_cluster_mean_psi.csv
   
4. Find single gene related to mean psi in a single cluster, do permutation to form a background distrbution,
   and substract genes variable across cell types.
   Input: SR_normed expression matrix from Seurat object.
   Output: correlated gene list
"""


############# 1. clustering to find local factors  #############
MN_se_psi_all= read.csv("/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/MN_psi/se/psi.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
NPC_se_psi_all= read.csv("/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/NPC_psi/se/psi.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
iPSC_se_psi_all= read.csv("/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/iPSC_psi/se/psi.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
MN_se_psi_all<-t(MN_se_psi_all)
NPC_se_psi_all<-t(NPC_se_psi_all)
iPSC_se_psi_all<-t(iPSC_se_psi_all)

#这里所用的函数在useful_functions.R
all_psi<-merge_by_rownames(MN_se_psi_all,NPC_se_psi_all)
all_psi<-merge_by_rownames(all_psi,iPSC_se_psi_all)

#filter116= read.csv("/Users/chenxinyi/Desktop/output/filtered116",sep=" ",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
#all_psi2<-all_psi[,which(colnames(all_psi) %in% filter116$x)]   #113个细胞，2104events
all_psi2<-all_psi[,!colnames(all_psi) %in% c("SRR4047260", "SRR4047278","SRR4047295","SRR4047315","SRR4047390","SRR4047332","SRR4047331", "SRR4047323")] 

# 初始有2104个AS
all_psi3<-all_psi2
all_psi3$nonna<-apply(all_psi2,1,sum_non_na)
hist(all_psi3$nonna,main="Histogram of nonNA cell for each AS",xlab="Number of nonNA cell for each AS")

all_psi_filt<-all_psi3[which(all_psi3$nonna>=80),]  #这样就有112个细胞和151个基因
all_psi_filt<-all_psi_filt[,-ncol(all_psi_filt)]
#剩下的NA用一个来自该AS的分布来补上
for (i in 1:length(all_psi_filt[,1])){
  len_null <- length(all_psi_filt[i,which(is.na(all_psi_filt[i,]))])
  all_psi_filt[i,which(is.na(all_psi_filt[i,]))] <- sample(as.numeric(all_psi_filt[i,which(!is.na(all_psi_filt[i,]))]))[1:len_null]
}

write.table (all_psi_filt, file = "/Users/chenxinyi/Desktop/output/all_se_psi_filt100cells.csv", row.names =TRUE, col.names =TRUE, quote =TRUE)

#另一个选项，直接计算平均欧式距离
AS_distance <- data.frame(matrix(NA,nrow(all_psi_filt),nrow(all_psi_filt)))
colnames(AS_distance)<-rownames(all_psi_filt)
rownames(AS_distance)<-rownames(all_psi_filt)

#一个计算距离和熵的宝藏555
#install.packages("philentropy")
library(philentropy)

for (row1 in 1:nrow(AS_distance)){
  for (col1 in 1:row1){
    AS_distance[row1,col1]<-calculate_dist_with_na(all_psi_filt[row1,],all_psi_filt[col1,])
  }
}

for (row1 in 1:nrow(AS_distance)){
  for (col1 in (row1+1):nrow(AS_distance)){
    AS_distance[row1,col1]<-AS_distance[col1,row1]
  }
}

write.table (AS_distance[1:537,1:537], file = "/Users/chenxinyi/Desktop/output/AS_distance.csv", row.names =TRUE, col.names =TRUE, quote =TRUE)


write.table (all_psi3, file = "/Users/chenxinyi/Desktop/all_se_psi_nonna.csv", row.names =TRUE, col.names =TRUE, quote =TRUE)

########## 2. find GO for different clusters ###########
label_pred<- read.csv("/Users/chenxinyi/Desktop/output/AS80_label.csv",row.names = 1,stringsAsFactors = FALSE) 
a<-data.frame(rownames(AS_distance),label_pred)
colnames(a)<-c("event_id","label_pred")
events$event_id<- gsub("[^0-9A-z]",".",events$event_id)
event_AS_distance<-merge(events,a,by="event_id",all=FALSE)
cluster0<-unique(event_AS_distance$isoform1_gene_name[which(event_AS_distance$label_pred==0)])
cluster1<-unique(event_AS_distance$isoform1_gene_name[which(event_AS_distance$label_pred==1)])
cluster2<-unique(event_AS_distance$isoform1_gene_name[which(event_AS_distance$label_pred==2)])
cluster3<-unique(event_AS_distance$isoform1_gene_name[which(event_AS_distance$label_pred==3)])

DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = cluster1,
                       keytype = "SYMBOL",
                       column = "ENTREZID")   
#有67个基因没有转换成功
DEG.entrez_id = na.omit(DEG.entrez_id)
#biological process 富集分析
#gene_mat2<-read.csv (file = "/Users/chenxinyi/Desktop/gene_mat2.csv",sep=" ",header = TRUE,row.names = 1,stringsAsFactors = FALSE)
ref.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = unique(event_AS_distance$isoform1_gene_name),
                       keytype = "SYMBOL",
                       column = "ENTREZID")   
#有779个基因没有转换成功(因为之前从ensemble转成基因名就已经失败了)
ref.entrez_id = na.omit(ref.entrez_id)

erich.go.BP_ref = enrichGO(gene = DEG.entrez_id,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID",
                           ont = "BP",
                           pvalueCutoff = 0.5,
                           qvalueCutoff = 0.5)#,
                           #universe = as.character(ref.entrez_id))  # universe: background genes

##分析完成后，作图
dotplot(erich.go.BP_ref)

#这里也可以用全部基因作为背景

#找到相应的 gene

all_se_psi<- read.csv("/Users/chenxinyi/Desktop/output/all_se_psi_nonna.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 


########## 3. use psi in a cluster of AS to identify relevant pathways ############
event_label= read.csv("/Users/chenxinyi/Desktop/output/event_label.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
length(event_label$event[which(event_label$label==0)])  #112个
length(event_label$event[which(event_label$label==1)])  #39个
grp1_list<-event_label$event[which(event_label$label==0)]
grp2_list<-event_label$event[which(event_label$label==1)]
grp1<-all_psi2[which(rownames(all_psi2) %in% grp1_list),]
grp2<-all_psi2[which(rownames(all_psi2) %in% grp2_list),]

#grp1<-t(grp1)
#grp1<-as.data.frame(grp1)
#grp1_orig<-grp1
#grp1$mean_psi<-0
#for (i in 1:length(grp1$mean_psi)){
#  grp1$mean_psi[i]<-mean(grp1_orig[i,])         #要计算new psi 吗？？？这里的NA好多
#}


########## 3.1 calculate global mean psi for selected genes (MN) ###############
#incorrect (look at recalculate_psi_less_coverage.R)
  
################ 3.2 NPC ##################

events= read.csv("/Users/chenxinyi/Desktop/output/NPC_mean_psi/events.csv",header = TRUE,stringsAsFactors = FALSE) 
junctions<-events[,c("event_id","junction13","junction12","junction23")]

junc_reads= read.csv("/Users/chenxinyi/Desktop/output/NPC_mean_psi/reads.csv",header = TRUE,stringsAsFactors = FALSE) 
head(junc_reads)

junctions$event_id<- gsub("[^0-9A-z]",".",junctions$event_id)

SRRlist<-unique(junc_reads$sample_id)
sum_junc<-data.frame(matrix(NA,3,length(SRRlist)))
colnames(sum_junc)<-SRRlist
rownames(sum_junc)<-c("junc13","junc12","junc23")
sum_junc_grp2<-sum_junc

#对于SRR1
for (j in 1:length(SRRlist)){
  junctions$read13<-0
  junctions$read12<-0
  junctions$read23<-0
  cell_junc<-junc_reads[which(junc_reads$sample_id==SRRlist[j]),]
  
  junction_clus1<-junctions[which(junctions$event_id %in% grp1_list),]   #提出了1169行
  for (i in 1:length(junction_clus1[,1])){   #不管正负链，永远是小的在前，大的在后的
    #junc_reads$sample_id[which(junc_reads$junction_id==junctions$junction13[i])]
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction13[i])])==1){
      if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction12[i])])==1){
        if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction23[i])])==1){
          junction_clus1$read13[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction13[i])]
          junction_clus1$read12[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction12[i])]
          junction_clus1$read23[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction23[i])]
        }
      } 
    }
  }
  adeq_junc<-which(junction_clus1$read13>=10 & junction_clus1$read12>=10 & junction_clus1$read23>=10)
  sum_junc[1,j]<-sum(junction_clus1[adeq_junc,"read13"])
  sum_junc[2,j]<-sum(junction_clus1[adeq_junc,"read12"])
  sum_junc[3,j]<-sum(junction_clus1[adeq_junc,"read23"])
  
  
  junction_clus2<-junctions[which(junctions$event_id %in% grp2_list),]   #提出了417行
  for (i in 1:length(junction_clus2[,1])){   #不管正负链，永远是小的在前，大的在后的
    #junc_reads$sample_id[which(junc_reads$junction_id==junctions$junction13[i])]
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction13[i])])==1){
      if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction12[i])])==1){
        if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction23[i])])==1){
          junction_clus2$read13[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction13[i])]
          junction_clus2$read12[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction12[i])]
          junction_clus2$read23[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction23[i])]
        }
      } 
    }
  }
  adeq_junc<-which(junction_clus2$read13>=10 & junction_clus2$read12>=10 & junction_clus2$read23>=10)
  sum_junc_grp2[1,j]<-sum(junction_clus2[adeq_junc,"read13"])
  sum_junc_grp2[2,j]<-sum(junction_clus2[adeq_junc,"read12"])
  sum_junc_grp2[3,j]<-sum(junction_clus2[adeq_junc,"read23"])
  
  print(j)
}

NPC_clus1_mean_psi<-as.data.frame(t(sum_junc))
NPC_clus1_mean_psi$mean_psi_clus1<-(NPC_clus1_mean_psi$junc12 + NPC_clus1_mean_psi$junc23)/(NPC_clus1_mean_psi$junc12 + NPC_clus1_mean_psi$junc23 + 2*NPC_clus1_mean_psi$junc13)
NPC_clus2_mean_psi<-as.data.frame(t(sum_junc_grp2))
NPC_clus2_mean_psi$mean_psi_clu2<-(NPC_clus2_mean_psi$junc12 + NPC_clus2_mean_psi$junc23)/(NPC_clus2_mean_psi$junc12 + NPC_clus2_mean_psi$junc23 + 2*NPC_clus2_mean_psi$junc13)

NPC_mean_psi_clustered<-merge(NPC_clus1_mean_psi,NPC_clus2_mean_psi,by=0,all=FALSE)
write.table (NPC_mean_psi_clustered, file = "/Users/chenxinyi/Desktop/NPC_cluster_mean_psi.csv", row.names =TRUE, col.names =TRUE, quote =TRUE)

#plot(NPC_mean_psi_clustered$mean_psi_clus1,NPC_mean_psi_clustered$mean_psi_clu2,xlab="mean psi for group1 genes (sum of junction reads)",ylab="mean psi for group2 genes (sum of junction reads)")

summary(lm(NPC_mean_psi_clustered$mean_psi_clu2 ~ NPC_mean_psi_clustered$mean_psi_clus1))

################ 3.3 iPSC #################

events= read.csv("/Users/chenxinyi/Desktop/output/iPSC_mean_psi/events.csv",header = TRUE,stringsAsFactors = FALSE) 
junctions<-events[,c("event_id","junction13","junction12","junction23")]

junc_reads= read.csv("/Users/chenxinyi/Desktop/output/iPSC_mean_psi/reads.csv",header = TRUE,stringsAsFactors = FALSE) 
head(junc_reads)

junctions$event_id<- gsub("[^0-9A-z]",".",junctions$event_id)

SRRlist<-unique(junc_reads$sample_id)
sum_junc<-data.frame(matrix(NA,3,length(SRRlist)))
colnames(sum_junc)<-SRRlist
rownames(sum_junc)<-c("junc13","junc12","junc23")
sum_junc_grp2<-sum_junc

#对于SRR1
for (j in 1:length(SRRlist)){
  junctions$read13<-0
  junctions$read12<-0
  junctions$read23<-0
  cell_junc<-junc_reads[which(junc_reads$sample_id==SRRlist[j]),]
  
  junction_clus1<-junctions[which(junctions$event_id %in% grp1_list),]   #提出了1169行
  for (i in 1:length(junction_clus1[,1])){   #不管正负链，永远是小的在前，大的在后的
    #junc_reads$sample_id[which(junc_reads$junction_id==junctions$junction13[i])]
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction13[i])])==1){
      if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction12[i])])==1){
        if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction23[i])])==1){
          junction_clus1$read13[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction13[i])]
          junction_clus1$read12[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction12[i])]
          junction_clus1$read23[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction23[i])]
        }
      } 
    }
  }
  adeq_junc<-which(junction_clus1$read13>=10 & junction_clus1$read12>=10 & junction_clus1$read23>=10)
  sum_junc[1,j]<-sum(junction_clus1[adeq_junc,"read13"])
  sum_junc[2,j]<-sum(junction_clus1[adeq_junc,"read12"])
  sum_junc[3,j]<-sum(junction_clus1[adeq_junc,"read23"])
  
  
  junction_clus2<-junctions[which(junctions$event_id %in% grp2_list),]   #提出了417行
  for (i in 1:length(junction_clus2[,1])){   #不管正负链，永远是小的在前，大的在后的
    #junc_reads$sample_id[which(junc_reads$junction_id==junctions$junction13[i])]
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction13[i])])==1){
      if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction12[i])])==1){
        if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction23[i])])==1){
          junction_clus2$read13[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction13[i])]
          junction_clus2$read12[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction12[i])]
          junction_clus2$read23[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction23[i])]
        }
      } 
    }
  }
  adeq_junc<-which(junction_clus2$read13>=10 & junction_clus2$read12>=10 & junction_clus2$read23>=10)
  sum_junc_grp2[1,j]<-sum(junction_clus2[adeq_junc,"read13"])
  sum_junc_grp2[2,j]<-sum(junction_clus2[adeq_junc,"read12"])
  sum_junc_grp2[3,j]<-sum(junction_clus2[adeq_junc,"read23"])
  
  print(j)
}

iPSC_clus1_mean_psi<-as.data.frame(t(sum_junc))
iPSC_clus1_mean_psi$mean_psi_clus1<-(iPSC_clus1_mean_psi$junc12 + iPSC_clus1_mean_psi$junc23)/(iPSC_clus1_mean_psi$junc12 + iPSC_clus1_mean_psi$junc23 + 2*iPSC_clus1_mean_psi$junc13)
iPSC_clus2_mean_psi<-as.data.frame(t(sum_junc_grp2))
iPSC_clus2_mean_psi$mean_psi_clu2<-(iPSC_clus2_mean_psi$junc12 + iPSC_clus2_mean_psi$junc23)/(iPSC_clus2_mean_psi$junc12 + iPSC_clus2_mean_psi$junc23 + 2*iPSC_clus2_mean_psi$junc13)

iPSC_mean_psi_clustered<-merge(iPSC_clus1_mean_psi,iPSC_clus2_mean_psi,by=0,all=FALSE)
write.table (iPSC_mean_psi_clustered, file = "/Users/chenxinyi/Desktop/iPSC_cluster_mean_psi.csv", row.names =TRUE, col.names =TRUE, quote =TRUE)

plot(iPSC_mean_psi_clustered$mean_psi_clus1,iPSC_mean_psi_clustered$mean_psi_clu2,xlab="mean psi for group1 genes (sum of junction reads)",ylab="mean psi for group2 genes (sum of junction reads)")

#summary(lm(NPC_mean_psi_clustered$mean_psi_clu2 ~ NPC_mean_psi_clustered$mean_psi_clus1))
MN_mean_psi_clustered$cell<-"MN"
NPC_mean_psi_clustered$cell<-"NPC"
iPSC_mean_psi_clustered$cell<-"iPSC"
a<-rbind(MN_mean_psi_clustered,NPC_mean_psi_clustered)
all_2clust<-rbind(a,iPSC_mean_psi_clustered)

library(ggplot2)
ggplot(all_2clust,aes(x=`mean_psi_clus1`,y=`mean_psi_clu2`,color=as.factor(cell))) + geom_point(shape=19) +
  geom_point() + 
  scale_colour_brewer(palette = "Set1") +
  #geom_smooth(method="lm",show.legend = TRUE,level=0.95) +
  #scale_fill_discrete(breaks=c("5_n1","5_n2","3_n1","3_n2")) +
  xlab("mean psi for cluster1 gene") + ylab("mean psi for cluster2 gene") 

############## 4. Try identify single genes related to these patterns #############
normed_expr<-as.data.frame(SR_normed@assays$RNA@data) #将一个稀疏矩阵转化成正常的矩阵
normed_expr<-t(normed_expr)
rownames(normed_expr)<-gsub("TPM_", "", rownames(normed_expr))
var_expr<-normed_expr[,which(colnames(normed_expr) %in% VariableFeatures(SR_normed))]   #取头2000个highly variable的基因去做分析

rownames(all_2clust)<-all_2clust$Row.names
comb<-merge(var_expr,all_2clust,by=0,all=FALSE)
rownames(comb)<-comb$Row.names
comb<-comb[,-1]

row<- data.frame(matrix(NA,2,length(colnames(comb))-10))
rownames(row)<-c("clust1","clust2")    # 都spearman correlation
colnames(row)<-colnames(comb)[1:(length(colnames(comb))-10)]

for (i in 1:(length(colnames(comb))-10)){
  row[1,i] <-cor(comb[,i],comb[,"mean_psi_clus1"],method = "spearman")
  row[2,i] <-cor(comb[which(!is.na(comb$mean_psi_clu2)),i],comb[which(!is.na(comb$mean_psi_clu2)),"mean_psi_clu2"],method = "spearman")
}
hist(as.numeric(row[1,]),breaks = 10,xlab = "spearman corr",main = "Histogram of spearman corr between cluster1 mean psi \nand gene expression")
hist(as.numeric(row[2,]),breaks = 10,xlab = "spearman corr",main = "Histogram of spearman corr between cluster2 mean psi \nand gene expression")

############# 4.2 Permutation ##################
iter<-500
row2<- data.frame(matrix(NA,iter,length(colnames(comb))-10))
#rownames(row2)<-"pearson_cor"
colnames(row2)<-colnames(comb)[1:(length(colnames(comb))-10)]
row3<-row2

for (j in 1:iter){
  #rand1<-sample(comb$mean_psi_clus1)
  rand2<-sample(comb[which(!is.na(comb$mean_psi_clu2)),"mean_psi_clu2"])
  for (i in 1:(length(colnames(comb))-10)){
    #row2[j,i] <-cor(comb[,i],rand1,method = "spearman")
    row3[j,i] <-cor(comb[which(!is.na(comb$mean_psi_clu2)),i],rand2,method = "spearman")
  }
}
library(reshape2)
agg2<-melt(row2)
agg3<-melt(row3)
hist(as.numeric(row[2,]),freq=FALSE,col=rgb(1,0,0,1/4),breaks = 20,xlim=range(-1,1),xlab = "spearman corr",add=T,main = "Histogram of spearman corr between randomized psi \nand gene expression")
hist(agg3$value,col=rgb(0,0,1,1/4),freq=FALSE,breaks = 20,xlim=range(-1,1),xlab = "spearman corr",main = "Histogram of spearman corr between cluter2 mean psi \nand gene expression")
legend("topleft", c("observed", "permutated"),fill=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)),bty="n")  #bty="n"是legend周围的边框不要

############ 4.3 把cell type之间的DEG过滤掉 ############
MN_markers <- FindMarkers(SR, ident.1 = "MN", ident.2 = "NPC", min.pct = 0.25)    #only.pos=FALSE by default
NPC_markers <- FindMarkers(SR, ident.1 = "NPC", ident.2 = "iPSC", min.pct = 0.25)    #only.pos=FALSE by default
iPSC_markers <- FindMarkers(SR, ident.1 = "iPSC", ident.2 = "MN", min.pct = 0.25)    #only.pos=FALSE by default
#特别多：11013，12096，12382
a<-union(rownames(MN_markers), rownames(NPC_markers))
DEG<- union(a,rownames(iPSC_markers))  #12603个DEG

############ 4.4 分析cluster1有什么显著的基因 ##############
#火山图画一下

vol_data<-data.frame(matrix(NA,(length(colnames(comb))-10),2))
rownames(vol_data)<-colnames(comb)[1:(length(colnames(comb))-10)]
colnames(vol_data)<-c("corr","pvalue")
for (i in 1:(length(colnames(comb))-10)){
  a<-cor.test(comb[,i],comb[,"mean_psi_clus1"],method = "pearson")
  vol_data[i,1] <-a$estimate
  vol_data[i,2] <-a$p.value
}
volcano<-subset(vol_data,select = c("corr","pvalue"))
threshold<-as.factor((abs(volcano$corr)>0.4) & volcano$pvalue<2.5e-5)   #其实最好是用bonferonni corrected P value
r03=ggplot(volcano,aes(corr,-log2(pvalue),colour=threshold))+geom_point()
r04=r03+labs(title="Volcanoplot (Bonferroni-corrected test)")+theme(plot.title = element_text(hjust = 0.5))+xlim(-1,1)
r05=r04+geom_vline(xintercept=c(-0.4,0.4),linetype="dotted",size=1)+geom_hline(yintercept=-log2(2.5e-5),col="blue")

library(ggrepel)    #可以让字与字自动错开
label<-rownames(volcano)
label[which((abs(volcano$corr)<=0.4) | volcano$pvalue>=2.5e-5)]<-NA
r06=r05+geom_text_repel(aes(corr, -log2(pvalue), label = label))

sign1<-rownames(volcano[which((abs(volcano$corr)>0.4) & volcano$pvalue<2.5e-5),])
setdiff(sign1,DEG) #属于显著基因，但又不是DEG的

############ 4.5 分析cluster2有什么显著的基因 ##############
vol_data<-data.frame(matrix(NA,(length(colnames(comb))-10),2))
rownames(vol_data)<-colnames(comb)[1:(length(colnames(comb))-10)]
colnames(vol_data)<-c("corr","pvalue")
for (i in 1:(length(colnames(comb))-10)){
  a<-cor.test(comb[which(!is.na(comb$mean_psi_clu2)),i],comb[which(!is.na(comb$mean_psi_clu2)),"mean_psi_clu2"],method = "pearson")
  vol_data[i,1] <-a$estimate
  vol_data[i,2] <-a$p.value
}
volcano<-subset(vol_data,select = c("corr","pvalue"))
threshold<-as.factor((abs(volcano$corr)>0.5) & volcano$pvalue<2.5e-5)   #其实最好是用bonferonni corrected P value
r03=ggplot(volcano,aes(corr,-log2(pvalue),colour=threshold))+geom_point()
r04=r03+labs(title="Volcanoplot for cluster2 (Bonferroni-corrected test)")+theme(plot.title = element_text(hjust = 0.5))+xlim(-1,1)
r05=r04+geom_vline(xintercept=c(-0.5,0.5),linetype="dotted",size=1)+geom_hline(yintercept=-log2(2.5e-5),col="blue")

library(ggrepel)    #可以让字与字自动错开
label<-rownames(volcano)
label[which((abs(volcano$corr)<=0.5) | volcano$pvalue>=2.5e-5)]<-NA
r06=r05+geom_text_repel(aes(corr, -log2(pvalue), label = label))
r06

sign2<-rownames(volcano[which((abs(volcano$corr)>0.5) & volcano$pvalue<2.5e-5),])
setdiff(sign2,DEG) #属于显著基因，但又不是DEG的

