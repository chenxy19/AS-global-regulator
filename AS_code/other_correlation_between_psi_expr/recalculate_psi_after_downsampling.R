"""
1. Calculate psi based on sum of junction reads
  In: events.csv, reads.csv from outrigger. The former is the mapping from alternative splicing
      events to junctions names, while the latter is the read count for all junctions.
  Out: (1) sum_junc.csv. the sum of junction13, junction12 and junction23 in all cells. 
       (2) junc_mat.csv. the count read for each junction (in the order of junction13, junction12 and junction23)
           call function:calculate_psi to convert into AS event list * sample list.

2. Downsampling: downsample raw reads/count to the sample with the least number of reads/counts
  and recalculate psi after downsampling.
  Downsampling avoids the technical effect introduced by uneven sequencing depth.

3. Calculate and plot two types of mean psi
  write into all_mean_psi_corrected0831.csv

4. Find the source of correlation betweem n_feature and gene-averaged psi.
  Draw the scatter plot between n_feature and gene-averaged psi (from all genes or the most highly expressed genes).
  Input: SR Seurat project.
  
5. find genes with high correlation with mean psi of two clusters of AS
  Input: grp1_list
  Output: cluster_mean_psi.csv
  
6. plot the number of reads in MN, NPC, iPSC
"""



########## 1. calculate psi with lower threshold (MN) ###############
for (prefix in c("mouseSC","mouseSC2")){
  event_file=paste("/Users/chenxinyi/Desktop/output/",prefix,"_events.csv",sep = "")
  events= read.csv(event_file,header = TRUE,stringsAsFactors = FALSE) 
  junctions<-events[,c("event_id","junction13","junction12","junction23")]    #其实unique只有8690,which is reasonable
  junctions<-unique(junctions)
  
  junc_reads_file=paste("/Users/chenxinyi/Desktop/output/",prefix,"_reads.csv",sep = "")
  junc_reads= read.csv(junc_reads_file,header = TRUE,stringsAsFactors = FALSE) 
  #head(junc_reads)
  
  junctions$event_id<- gsub("[^0-9A-z]",".",junctions$event_id)
  
  SRRlist<-unique(junc_reads$sample_id)
  sum_junc<-data.frame(matrix(NA,3,length(SRRlist)))
  colnames(sum_junc)<-SRRlist
  rownames(sum_junc)<-c("junc13","junc12","junc23")
  ASlist<-junctions$event_id
  
  #junc_mat<-data.frame(matrix(0,3*length(ASlist),length(SRRlist)))
  #colnames(junc_mat)<-SRRlist
  #for (i in 1:length(ASlist)){
  #  rownames(junc_mat)[3*i-2] <- paste(ASlist[i],"junc13",sep="_")
  #  rownames(junc_mat)[3*i-1] <- paste(ASlist[i],"junc12",sep="_")
  #  rownames(junc_mat)[3*i] <- paste(ASlist[i],"junc23",sep="_")
  #}
  
  #对于SRR1
  for (j in 1:length(SRRlist)){
    cell_junc<-junc_reads[which(junc_reads$sample_id==SRRlist[j]),]
    junctions$read13<-0
    junctions$read12<-0
    junctions$read23<-0
    for (i in 1:length(junctions[,1])){   #不管正负链，永远是小的在前，大的在后的
      #junc_reads$sample_id[which(junc_reads$junction_id==junctions$junction13[i])]
      #不需要每个junction必须有reads!这里太苛刻了！
      tp<-cell_junc$reads[which(cell_junc$junction_id==junctions$junction13[i])]
      if (length(tp)==1){
        junctions$read13[i]<-tp
      }
      if (length(tp)>1){
        print(paste(junctions$junction13[i]," multiple occur!",sep=""))
      }
      
      tp<-cell_junc$reads[which(cell_junc$junction_id==junctions$junction12[i])]
      if (length(tp)==1){
        junctions$read12[i]<-tp
      }
      
      tp<-cell_junc$reads[which(cell_junc$junction_id==junctions$junction23[i])]
      if (length(tp)==1){
        junctions$read23[i]<-tp
      }
  }
  
    adeq_junc<-which(junctions$read13 + junctions$read12 + junctions$read23 >= 10)  #从134个升到169，提升并不显著
    sum_junc[1,j]<-sum(junctions[adeq_junc,"read13"])
    sum_junc[2,j]<-sum(junctions[adeq_junc,"read12"])
    sum_junc[3,j]<-sum(junctions[adeq_junc,"read23"])
    #for (k in adeq_junc){
    #  MN_psi[k,SRRlist[j]]<-(junctions$read12[k] + junctions$read23[k])/(junctions$read12[k] + junctions$read23[k] + 2*junctions$read13[k])
    #}
    print(SRRlist[j])
    #head(junc_mat[,j])
  }
  
  sum_junc_file=paste("/Users/chenxinyi/Desktop/output/",prefix,"_sum_junc.csv",sep="")
  write.table (sum_junc, file = sum_junc_file, row.names =TRUE, col.names =TRUE, quote =TRUE)

  #MN_globl_mean_psi<-as.data.frame(t(sum_junc))
  #MN_globl_mean_psi$mean_psi<-(MN_globl_mean_psi$junc12 + MN_globl_mean_psi$junc23)/(MN_globl_mean_psi$junc12 + MN_globl_mean_psi$junc23 + 2*MN_globl_mean_psi$junc13)
}

#发现重复的isoform好多，决定前面event就可以unique一下

########## 2. downsampling ################
"""
calculate_psi: a function in useful_functions.R. 
  Input: junc_mat.csv.
  Output: a dataframe of AS list * sample list.
"""
iPSC_junc_mat<-junc_mat
iPSC_psi_junc_cell <- calculate_psi("iPSC",junc_mat)

junc_mat_file=paste("/Users/chenxinyi/Desktop/output/MN_junc_mat.csv",sep="")
MN_junc_mat= read.csv(junc_mat_file,sep=" ",header = TRUE,stringsAsFactors = FALSE) 
MN_psi_junc_cell <- calculate_psi("MN",MN_junc_mat)

junc_mat_file=paste("/Users/chenxinyi/Desktop/output/NPC_junc_mat.csv",sep="")
NPC_junc_mat= read.csv(junc_mat_file,sep=" ",header = TRUE,stringsAsFactors = FALSE) 
NPC_psi_junc_cell <- calculate_psi("NPC",NPC_junc_mat)

#应该只留过了QC的116个
filtered <- read.csv("/Users/chenxinyi/Desktop/output/filtered116",sep=" ", header = TRUE,stringsAsFactors = FALSE) 

min_reads<-10000000
for (i in which(colnames(iPSC_junc_mat) %in% filtered$x)){
  if (min_reads > sum(iPSC_junc_mat[,i])){
    min_reads<-sum(iPSC_junc_mat[,i])
  }
}
# iPSC 的最小是 1750694
for (i in which(colnames(MN_junc_mat) %in% filtered$x)){
  if (min_reads > sum(MN_junc_mat[,i])){
    min_reads<-sum(MN_junc_mat[,i])
  }
}
#MN 最小是 165703
for (i in which(colnames(NPC_junc_mat) %in% filtered$x)){
  if (min_reads > sum(NPC_junc_mat[,i])){
    min_reads<-sum(NPC_junc_mat[,i])
  }
}
#NPC最小是506654

#都downsample到165703

#source("https://bioconductor.org/biocLite.R")
#biocLite("DropletUtils")
library(DropletUtils)
# downsample is a function in useful_functions.R
MN_junc_mat_QC_downsampled <- downsample(MN_junc_mat)
NPC_junc_mat_QC_downsampled <- downsample(NPC_junc_mat)
iPSC_junc_mat_QC_downsampled <- downsample(iPSC_junc_mat)

MN_psi_downsampled<-calculate_psi("MN",MN_junc_mat_QC_downsampled,SRRlist=intersect(colnames(MN_junc_mat),filtered$x))
NPC_psi_downsampled<-calculate_psi("NPC",NPC_junc_mat_QC_downsampled,SRRlist=intersect(colnames(NPC_junc_mat),filtered$x))
iPSC_psi_downsampled<-calculate_psi("iPSC",iPSC_junc_mat_QC_downsampled,SRRlist=intersect(colnames(iPSC_junc_mat),filtered$x))


all_psi_downsampled<-merge_by_rownames(MN_psi_downsampled,NPC_psi_downsampled)
all_psi_downsampled<-merge_by_rownames(all_psi_downsampled,iPSC_psi_downsampled)
#113 cell* 2104 events,超过一般都是NA。。。


num_nonna<-apply(all_psi_downsampled,1,sum_non_na)  #第一个是矩阵，第二个是对行操作（2 是对列操作），第三个是函数
hist(num_nonna,xlab="number of non NA cells",main="Histogram of the number of cells\n an AS is detected")

all_psi_downsampled_high <- all_psi_downsampled[which(num_nonna>110),]
for (i in 1:nrow(all_psi_downsampled_high)){
  if(sum_non_na(all_psi_downsampled_high[i,])!=ncol(all_psi_downsampled_high)){
    all_psi_downsampled_high[i,which(is.na(all_psi_downsampled_high[i,]))] <- mean(as.numeric(all_psi_downsampled_high[i,which(!is.na(all_psi_downsampled_high[i,]))]))
  }
}

#怎么衡量一个cell population的异质性？
#选取全部检测到的;有74个基因在超过110个细胞中被检测到，用于tSNE,剩余部分用该AS的平均值补上
#tSNE
library(Rtsne)  #perplexity是计算时考虑多少个neighbor，不应该超过总样本数
tsne<- Rtsne(t(all_psi_downsampled_high), dims = 2, perplexity=10, verbose=TRUE, max_iter = 1000)
Labels<-rep(c("MN","NPC","iPSC"),times=c(ncol(MN_psi_downsampled),ncol(NPC_psi_downsampled),ncol(iPSC_psi_downsampled)))
Labels<-as.factor(Labels)
colors = rainbow(length(unique(Labels)))  #每个点按照细胞类型上色

plot(tsne$Y, xlab="tsne_1",ylab="tsne_2",main="tsne after downsampling to 165703 \njunction reads/cell (perplexity=10)",col=colors[Labels])
legend("bottomright",levels(Labels),col  = colors,pch=15)
#text(tsne$Y, labels=Labels, col=colors[Labels])  #每个点都用plot来表示


#对downsampling之前的做一下
all_psi_junc_cell<-merge_by_rownames(MN_psi_junc_cell,NPC_psi_junc_cell)
all_psi_junc_cell<-merge_by_rownames(all_psi_junc_cell,iPSC_psi_junc_cell)
all_psi_junc_QC<-all_psi_junc_cell[,which(colnames(all_psi_junc_cell) %in% filtered$x)]
num_nonna_original<-apply(all_psi_junc_QC,1,sum_non_na)  #第一个是矩阵，第二个是对行操作（2 是对列操作），第三个是函数
hist(num_nonna_original,xlab="number of non NA cells",main="Histogram of the number of cells\n an AS is detected before downsampling")

all_psi_junc_QC_high <- all_psi_junc_QC[which(num_nonna_original>110),]  #  这个可以用前113维来做tSNE,包含的信息更多一些
for (i in 1:nrow(all_psi_junc_QC_high)){
  if(sum_non_na(all_psi_junc_QC_high[i,])!=ncol(all_psi_junc_QC_high)){
    all_psi_junc_QC_high[i,which(is.na(all_psi_junc_QC_high[i,]))] <- mean(as.numeric(all_psi_junc_QC_high[i,which(!is.na(all_psi_junc_QC_high[i,]))]))
  }
}

tsne<- Rtsne(t(all_psi_junc_QC_high), dims = 2, perplexity=10, verbose=TRUE, max_iter = 3000)
Labels<-rep(c("MN","NPC","iPSC"),times=c(ncol(MN_psi_downsampled),ncol(NPC_psi_downsampled),ncol(iPSC_psi_downsampled)))
Labels<-as.factor(Labels)
colors = rainbow(length(unique(Labels)))  #每个点按照细胞类型上色

plot(tsne$Y, xlab="tsne_1",ylab="tsne_2",main="tsne before downsampling (perplexity=10)",col=colors[Labels])
legend("bottomright",levels(Labels),col  = colors,pch=15)


############ 3. plot old psi with new psi ##############
"""
calculate_gene_average_psi: a function in useful_functions.R. 
  Input: se_bi_psi.csv. sample * AS events, show psi for all binary AS events.
  Output: a new matrix with the gene-averaged psi for all samples.
  
merge_into_all_matrix: a function in useful_functions.R. 
  merging average psi with the original psi matrix.
"""
MN_se_psi_mean<-calculate_gene_average_psi("MN")
NPC_se_psi_mean<-calculate_gene_average_psi("NPC")
iPSC_se_psi_mean<-calculate_gene_average_psi("iPSC")

colnames(all_globl_mean_psi)[4]<-"glbl_mean_psi"

MN_globl_mean_psi<-merge_into_all_matrix("MN",MN_se_psi_mean)
NPC_globl_mean_psi<-merge_into_all_matrix("NPC",NPC_se_psi_mean)
iPSC_globl_mean_psi<-merge_into_all_matrix("iPSC",iPSC_se_psi_mean)
all_psi_new<-rbind(MN_globl_mean_psi,NPC_globl_mean_psi,iPSC_globl_mean_psi)   #三个东西同时rbind也是可以的！
write.table (all_psi_new, file = "/Users/chenxinyi/Desktop/output/all_mean_psi_corrected0831.csv", row.names =TRUE, col.names =TRUE, quote =TRUE)

ggplot(all_psi_new,aes(x=`glbl_mean_psi`,y=`mean_psi_genes`,color=as.factor(cell_type))) + geom_point(shape=19) +
  geom_point() + 
  scale_colour_brewer(palette = "Set1") +
  geom_smooth(method="lm",show.legend = TRUE,level=0.95) +
  #scale_fill_discrete(breaks=c("5_n1","5_n2","3_n1","3_n2")) +
  xlab("Global mean psi (by addition of junction reads)") + ylab("Average of psi for genes") 

cor.test(all_psi_new[which(all_psi_new$cell_type=="iPSC"),"mean_psi_genes"], all_psi_new[which(all_psi_new$cell_type=="iPSC"),"glbl_mean_psi"],method="spearman")

########### 4. 看n_feature和new psi有什么关系（画在一张图上）##########
iPSC_se_psi= read.csv("/Users/chenxinyi/Desktop/output/iPSC_se_bi_psi.csv",header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
       #295和315是pool，295 QC时没给去掉，说明QC时总counts数还是要看的
iPSC_se_psi<-iPSC_se_psi[!rownames(iPSC_se_psi) %in% c("SRR4047260", "SRR4047278","SRR4047295","SRR4047315","SRR4047390","SRR4047332","SRR4047331", "SRR4047323"),] #NPC 31

MN_se_psi_mean<-MN_se_psi
MN_se_psi_mean$mean_psi<-0
MN_se_psi_mean$partial_mean<-0

cell_number_each_AS<-apply(MN_se_psi,2,sum_non_na)   #1 是对行进行操作
hist(as.numeric(cell_number_each_AS),main="number of cells each AS is detected in",xlab = "number of cells" )

for (i in 1:length(MN_se_psi[,1])){
  MN_se_psi_mean$mean_psi[i]<-mean(as.numeric(MN_se_psi[i,which(!is.na(MN_se_psi[i,]))]))   #na.rm=TRUE
}

nfeature<-as.data.frame(SR$nFeature_RNA)
rownames(nfeature)<-gsub("TPM_", "", rownames(nfeature))

layout(matrix(1:6, 2, 3,byrow=T))

for (j in 3:8){
  which_gene_exceed <- which(as.numeric(cell_number_each_AS) > 5*j)  #top 233/192
  MN_se_psi_partial <- MN_se_psi[,which_gene_exceed]
  for (i in 1:length(MN_se_psi[,1])){
    MN_se_psi_mean$partial_mean[i]<-mean(as.numeric(MN_se_psi_partial[i,which(!is.na(MN_se_psi_partial[i,])) ]))
  }
  MN_psi_nfeature<-merge_by_rownames(MN_se_psi_mean, nfeature)
  x<-MN_psi_nfeature$`SR$nFeature_RNA`
  y<-MN_psi_nfeature$partial_mean
  plot(x,y,main=paste("Pearson corr=",round(cor.test(x,y)$estimate,5),"\np-value=",round(cor.test(x,y)$p.value,5),sep=""),
       xlab="n_feature",ylab=paste("mean_psi_above",as.character(5*j),"cells",sep=" "))
  LM<-lm(y~x)
  abline(LM)
  summary(LM)
}

layout(matrix(1:1, 2, 3,byrow=T))

x<-MN_psi_nfeature$`SR$nFeature_RNA`
y<-MN_psi_nfeature$mean_psi
plot(x,y,main=paste("Pearson corr=",round(cor.test(x,y)$estimate,5),"\np-value=",round(cor.test(x,y)$p.value,5),sep=""),
     xlab="n_feature",ylab="mean_psi_all")
LM<-lm(y~x)
abline(LM)
summary(LM)

all_psi_nfeature<-merge(all_psi_new,nfeature,by=0,all=FALSE)
rownames(all_psi_nfeature)<-all_psi_nfeature$Row.names
all_psi_nfeature<-all_psi_nfeature[,-1]

ggplot(all_psi_nfeature,aes(x=`SR$nFeature_RNA`,y=`mean_psi_genes`,color=as.factor(cell_type))) + geom_point(shape=19) +
  geom_point() + 
  scale_colour_brewer(palette = "Set1") +
  geom_smooth(method="lm",show.legend = TRUE,level=0.95) +
  #scale_fill_discrete(breaks=c("5_n1","5_n2","3_n1","3_n2")) +
  xlab("Number of features") + ylab("mean psi (by averaging psi for each gene)") 

cor.test(all_psi_nfeature[which(all_psi_nfeature$cell_type=="MN"),"SR$nFeature_RNA"], all_psi_nfeature[which(all_psi_nfeature$cell_type=="MN"),"mean_psi_genes"],method="spearman")

################ 5. find genes with high correlation #############
#取出4000个high variable gene，把DEG过滤掉
SR_normed <- FindVariableFeatures(SR_normed, selection.method = "vst", nfeatures = 4000)
VariableFeaturePlot(SR_normed)
LabelPoints(plot = plot1, points = top10,repel = TRUE,xnudge =0, ynudge=0)  #repel是为了字和字之间的间隔
var_expr<-normed_expr[,which(colnames(normed_expr) %in% VariableFeatures(SR_normed))]   #取头2000个highly variable的基因去做分析
var_expr<-var_expr[,which(colnames(var_expr) %in% setdiff(colnames(normed_expr),DEG))]  # 只剩605个了

comb<-merge(var_expr,all_psi_new,by=0,all=FALSE)
rownames(comb)<-comb$Row.names
comb<-comb[,-1]


########## 5.1 按2个cluster 来计算 ###########
events= read.csv("/Users/chenxinyi/Desktop/output/MN_mean_psi/events.csv",header = TRUE,stringsAsFactors = FALSE) 
junctions<-events[,c("event_id","junction13","junction12","junction23")]    #其实unique只有8690,which is reasonable
junctions<-unique(junctions)

junc_reads= read.csv("/Users/chenxinyi/Desktop/output/MN_mean_psi/reads.csv",header = TRUE,stringsAsFactors = FALSE) 
head(junc_reads)

junctions$event_id<- gsub("[^0-9A-z]",".",junctions$event_id)

SRRlist<-unique(junc_reads$sample_id)
sum_junc<-data.frame(matrix(NA,3,length(SRRlist)))
colnames(sum_junc)<-SRRlist
rownames(sum_junc)<-c("junc13","junc12","junc23")
ASlist<-junctions$event_id
sum_junc_grp2<-sum_junc

#MN_psi<-data.frame(matrix(NA,length(ASlist),length(SRRlist)))
#colnames(MN_psi)<-SRRlist
#MN_psi<-cbind(junctions,MN_psi)

#对于SRR1
for (j in 1:length(SRRlist)){
  junctions$read13<-0
  junctions$read12<-0
  junctions$read23<-0
  cell_junc<-junc_reads[which(junc_reads$sample_id==SRRlist[j]),]
  
  junction_clus1<-junctions[which(junctions$event_id %in% grp1_list),]   #提出了1169行
  for (i in 1:length(junction_clus1[,1])){   #不管正负链，永远是小的在前，大的在后的
    #junc_reads$sample_id[which(junc_reads$junction_id==junctions$junction13[i])]
    #不需要每个junction必须有reads!这里太苛刻了！
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction13[i])])==1){
      junction_clus1$read13[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction13[i])]
    }
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction12[i])])==1){
      junction_clus1$read12[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction12[i])]
    }
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction23[i])])==1){
      junction_clus1$read23[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction23[i])]
    }
  }
  
  adeq_junc<-which(junction_clus1$read13 + junction_clus1$read12 + junction_clus1$read23 >= 10)  #从134个升到169，提升并不显著
  sum_junc[1,j]<-sum(junction_clus1[adeq_junc,"read13"])
  sum_junc[2,j]<-sum(junction_clus1[adeq_junc,"read12"])
  sum_junc[3,j]<-sum(junction_clus1[adeq_junc,"read23"])
  #for (k in adeq_junc){
  #  MN_psi[k,SRRlist[j]]<-(junctions$read12[k] + junctions$read23[k])/(junctions$read12[k] + junctions$read23[k] + 2*junctions$read13[k])
  #}
  
  junction_clus2<-junctions[which(junctions$event_id %in% grp2_list),]   #提出了1169行
  for (i in 1:length(junction_clus2[,1])){   #不管正负链，永远是小的在前，大的在后的
    #junc_reads$sample_id[which(junc_reads$junction_id==junctions$junction13[i])]
    #不需要每个junction必须有reads!这里太苛刻了！
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction13[i])])==1){
      junction_clus2$read13[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction13[i])]
    }
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction12[i])])==1){
      junction_clus2$read12[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction12[i])]
    }
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction23[i])])==1){
      junction_clus2$read23[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction23[i])]
    }
  }
  
  adeq_junc<-which(junction_clus2$read13 + junction_clus2$read12 + junction_clus2$read23 >= 10)  #从134个升到169，提升并不显著
  sum_junc_grp2[1,j]<-sum(junction_clus2[adeq_junc,"read13"])
  sum_junc_grp2[2,j]<-sum(junction_clus2[adeq_junc,"read12"])
  sum_junc_grp2[3,j]<-sum(junction_clus2[adeq_junc,"read23"])
  
  
  print(j)
}

MN_clus1_mean_psi<-as.data.frame(t(sum_junc))
MN_clus1_mean_psi$mean_psi_clus1<-(MN_clus1_mean_psi$junc12 + MN_clus1_mean_psi$junc23)/(MN_clus1_mean_psi$junc12 + MN_clus1_mean_psi$junc23 + 2*MN_clus1_mean_psi$junc13)
MN_clus2_mean_psi<-as.data.frame(t(sum_junc_grp2))
MN_clus2_mean_psi$mean_psi_clus2<-(MN_clus2_mean_psi$junc12 + MN_clus2_mean_psi$junc23)/(MN_clus2_mean_psi$junc12 + MN_clus2_mean_psi$junc23 + 2*MN_clus2_mean_psi$junc13)

MN_mean_psi_clustered<-merge(MN_clus1_mean_psi,MN_clus2_mean_psi,by=0,all=FALSE)
write.table (MN_mean_psi_clustered, file = "/Users/chenxinyi/Desktop/output/MN_cluster_mean_psi.csv", row.names =TRUE, col.names =TRUE, quote =TRUE)

plot(MN_mean_psi_clustered$mean_psi_clus1,MN_mean_psi_clustered$mean_psi_clus2)

########## 5.2 按2个cluster 来计算 NPC ###########
events= read.csv("/Users/chenxinyi/Desktop/output/NPC_mean_psi/events.csv",header = TRUE,stringsAsFactors = FALSE) 
junctions<-events[,c("event_id","junction13","junction12","junction23")]    #其实unique只有8690,which is reasonable
junctions<-unique(junctions)

junc_reads= read.csv("/Users/chenxinyi/Desktop/output/NPC_mean_psi/reads.csv",header = TRUE,stringsAsFactors = FALSE) 
head(junc_reads)

junctions$event_id<- gsub("[^0-9A-z]",".",junctions$event_id)

SRRlist<-unique(junc_reads$sample_id)
sum_junc<-data.frame(matrix(NA,3,length(SRRlist)))
colnames(sum_junc)<-SRRlist
rownames(sum_junc)<-c("junc13","junc12","junc23")
ASlist<-junctions$event_id
sum_junc_grp2<-sum_junc

#MN_psi<-data.frame(matrix(NA,length(ASlist),length(SRRlist)))
#colnames(MN_psi)<-SRRlist
#MN_psi<-cbind(junctions,MN_psi)

#对于SRR1
for (j in 1:length(SRRlist)){
  junctions$read13<-0
  junctions$read12<-0
  junctions$read23<-0
  cell_junc<-junc_reads[which(junc_reads$sample_id==SRRlist[j]),]
  
  junction_clus1<-junctions[which(junctions$event_id %in% grp1_list),]   #提出了1169行
  for (i in 1:length(junction_clus1[,1])){   #不管正负链，永远是小的在前，大的在后的
    #junc_reads$sample_id[which(junc_reads$junction_id==junctions$junction13[i])]
    #不需要每个junction必须有reads!这里太苛刻了！
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction13[i])])==1){
      junction_clus1$read13[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction13[i])]
    }
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction12[i])])==1){
      junction_clus1$read12[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction12[i])]
    }
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction23[i])])==1){
      junction_clus1$read23[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction23[i])]
    }
  }
  
  adeq_junc<-which(junction_clus1$read13 + junction_clus1$read12 + junction_clus1$read23 >= 10)  #从134个升到169，提升并不显著
  sum_junc[1,j]<-sum(junction_clus1[adeq_junc,"read13"])
  sum_junc[2,j]<-sum(junction_clus1[adeq_junc,"read12"])
  sum_junc[3,j]<-sum(junction_clus1[adeq_junc,"read23"])
  #for (k in adeq_junc){
  #  MN_psi[k,SRRlist[j]]<-(junctions$read12[k] + junctions$read23[k])/(junctions$read12[k] + junctions$read23[k] + 2*junctions$read13[k])
  #}
  
  junction_clus2<-junctions[which(junctions$event_id %in% grp2_list),]   #提出了1169行
  for (i in 1:length(junction_clus2[,1])){   #不管正负链，永远是小的在前，大的在后的
    #junc_reads$sample_id[which(junc_reads$junction_id==junctions$junction13[i])]
    #不需要每个junction必须有reads!这里太苛刻了！
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction13[i])])==1){
      junction_clus2$read13[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction13[i])]
    }
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction12[i])])==1){
      junction_clus2$read12[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction12[i])]
    }
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction23[i])])==1){
      junction_clus2$read23[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction23[i])]
    }
  }
  
  adeq_junc<-which(junction_clus2$read13 + junction_clus2$read12 + junction_clus2$read23 >= 10)  #从134个升到169，提升并不显著
  sum_junc_grp2[1,j]<-sum(junction_clus2[adeq_junc,"read13"])
  sum_junc_grp2[2,j]<-sum(junction_clus2[adeq_junc,"read12"])
  sum_junc_grp2[3,j]<-sum(junction_clus2[adeq_junc,"read23"])
  
  
  print(j)
}

NPC_clus1_mean_psi<-as.data.frame(t(sum_junc))
NPC_clus1_mean_psi$mean_psi_clus1<-(NPC_clus1_mean_psi$junc12 + NPC_clus1_mean_psi$junc23)/(NPC_clus1_mean_psi$junc12 + NPC_clus1_mean_psi$junc23 + 2*NPC_clus1_mean_psi$junc13)
NPC_clus2_mean_psi<-as.data.frame(t(sum_junc_grp2))
NPC_clus2_mean_psi$mean_psi_clus2<-(NPC_clus2_mean_psi$junc12 + NPC_clus2_mean_psi$junc23)/(NPC_clus2_mean_psi$junc12 + NPC_clus2_mean_psi$junc23 + 2*NPC_clus2_mean_psi$junc13)

NPC_mean_psi_clustered<-merge(NPC_clus1_mean_psi,NPC_clus2_mean_psi,by=0,all=FALSE)
write.table (NPC_mean_psi_clustered, file = "/Users/chenxinyi/Desktop/output/NPC_cluster_mean_psi.csv", row.names =TRUE, col.names =TRUE, quote =TRUE)

plot(NPC_mean_psi_clustered$mean_psi_clus1,NPC_mean_psi_clustered$mean_psi_clus2)

########## 5.3 按2个cluster 来计算 iPSC ###########
events= read.csv("/Users/chenxinyi/Desktop/output/iPSC_mean_psi/events.csv",header = TRUE,stringsAsFactors = FALSE) 
junctions<-events[,c("event_id","junction13","junction12","junction23")]    #其实unique只有8690,which is reasonable
junctions<-unique(junctions)

junc_reads= read.csv("/Users/chenxinyi/Desktop/output/iPSC_mean_psi/reads.csv",header = TRUE,stringsAsFactors = FALSE) 
head(junc_reads)

junctions$event_id<- gsub("[^0-9A-z]",".",junctions$event_id)

SRRlist<-unique(junc_reads$sample_id)
sum_junc<-data.frame(matrix(NA,3,length(SRRlist)))
colnames(sum_junc)<-SRRlist
rownames(sum_junc)<-c("junc13","junc12","junc23")
ASlist<-junctions$event_id
sum_junc_grp2<-sum_junc

#MN_psi<-data.frame(matrix(NA,length(ASlist),length(SRRlist)))
#colnames(MN_psi)<-SRRlist
#MN_psi<-cbind(junctions,MN_psi)

#对于SRR1
for (j in 1:length(SRRlist)){
  junctions$read13<-0
  junctions$read12<-0
  junctions$read23<-0
  cell_junc<-junc_reads[which(junc_reads$sample_id==SRRlist[j]),]
  
  junction_clus1<-junctions[which(junctions$event_id %in% grp1_list),]   #提出了1169行
  for (i in 1:length(junction_clus1[,1])){   #不管正负链，永远是小的在前，大的在后的
    #junc_reads$sample_id[which(junc_reads$junction_id==junctions$junction13[i])]
    #不需要每个junction必须有reads!这里太苛刻了！
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction13[i])])==1){
      junction_clus1$read13[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction13[i])]
    }
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction12[i])])==1){
      junction_clus1$read12[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction12[i])]
    }
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction23[i])])==1){
      junction_clus1$read23[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus1$junction23[i])]
    }
  }
  
  adeq_junc<-which(junction_clus1$read13 + junction_clus1$read12 + junction_clus1$read23 >= 10)  #从134个升到169，提升并不显著
  sum_junc[1,j]<-sum(junction_clus1[adeq_junc,"read13"])
  sum_junc[2,j]<-sum(junction_clus1[adeq_junc,"read12"])
  sum_junc[3,j]<-sum(junction_clus1[adeq_junc,"read23"])
  #for (k in adeq_junc){
  #  MN_psi[k,SRRlist[j]]<-(junctions$read12[k] + junctions$read23[k])/(junctions$read12[k] + junctions$read23[k] + 2*junctions$read13[k])
  #}
  
  junction_clus2<-junctions[which(junctions$event_id %in% grp2_list),]   #提出了1169行
  for (i in 1:length(junction_clus2[,1])){   #不管正负链，永远是小的在前，大的在后的
    #junc_reads$sample_id[which(junc_reads$junction_id==junctions$junction13[i])]
    #不需要每个junction必须有reads!这里太苛刻了！
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction13[i])])==1){
      junction_clus2$read13[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction13[i])]
    }
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction12[i])])==1){
      junction_clus2$read12[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction12[i])]
    }
    if (length(cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction23[i])])==1){
      junction_clus2$read23[i]<-cell_junc$reads[which(cell_junc$junction_id==junction_clus2$junction23[i])]
    }
  }
  
  adeq_junc<-which(junction_clus2$read13 + junction_clus2$read12 + junction_clus2$read23 >= 10)  #从134个升到169，提升并不显著
  sum_junc_grp2[1,j]<-sum(junction_clus2[adeq_junc,"read13"])
  sum_junc_grp2[2,j]<-sum(junction_clus2[adeq_junc,"read12"])
  sum_junc_grp2[3,j]<-sum(junction_clus2[adeq_junc,"read23"])
  
  
  print(j)
}

iPSC_clus1_mean_psi<-as.data.frame(t(sum_junc))
iPSC_clus1_mean_psi$mean_psi_clus1<-(iPSC_clus1_mean_psi$junc12 + iPSC_clus1_mean_psi$junc23)/(iPSC_clus1_mean_psi$junc12 + iPSC_clus1_mean_psi$junc23 + 2*iPSC_clus1_mean_psi$junc13)
iPSC_clus2_mean_psi<-as.data.frame(t(sum_junc_grp2))
iPSC_clus2_mean_psi$mean_psi_clus2<-(iPSC_clus2_mean_psi$junc12 + iPSC_clus2_mean_psi$junc23)/(iPSC_clus2_mean_psi$junc12 + iPSC_clus2_mean_psi$junc23 + 2*iPSC_clus2_mean_psi$junc13)

iPSC_mean_psi_clustered<-merge(iPSC_clus1_mean_psi,iPSC_clus2_mean_psi,by=0,all=FALSE)
write.table (iPSC_mean_psi_clustered, file = "/Users/chenxinyi/Desktop/output/iPSC_cluster_mean_psi.csv", row.names =TRUE, col.names =TRUE, quote =TRUE)

plot(iPSC_mean_psi_clustered$mean_psi_clus1,iPSC_mean_psi_clustered$mean_psi_clus2)

MN_mean_psi_clustered$cell_type<-"MN"
NPC_mean_psi_clustered$cell_type<-"NPC"
iPSC_mean_psi_clustered$cell_type<-"iPSC"
rownames(MN_mean_psi_clustered)<-MN_mean_psi_clustered$Row.names
rownames(NPC_mean_psi_clustered)<-NPC_mean_psi_clustered$Row.names
rownames(iPSC_mean_psi_clustered)<-iPSC_mean_psi_clustered$Row.names
all_cluster_psi<-rbind(MN_mean_psi_clustered,NPC_mean_psi_clustered,iPSC_mean_psi_clustered)
all_cluster_psi<-all_cluster_psi[,c("mean_psi_clus1","mean_psi_clus2","cell_type")]

ggplot(all_cluster_psi,aes(x=`mean_psi_clus1`,y=`mean_psi_clus2`,color=as.factor(cell_type))) + geom_point(shape=19) +
  geom_point() + 
  scale_colour_brewer(palette = "Set1") +
  geom_smooth(method="lm",show.legend = TRUE,level=0.95) +
  #scale_fill_discrete(breaks=c("5_n1","5_n2","3_n1","3_n2")) +
  xlab("cluster1 AS mean psi (by addition of junction reads)") + ylab("cluster2 AS mean psi (by addition of junction reads)") 

cor.test(all_cluster_psi[which(all_cluster_psi$cell_type=="iPSC"),"mean_psi_clus1"], all_cluster_psi[which(all_cluster_psi$cell_type=="iPSC"),"mean_psi_clus2"],method="spearman")


######### 5.4 FInd genes correlated with cluster1 AS ##########
comb<-merge(var_expr,all_cluster_psi,by=0,all=FALSE)
rownames(comb)<-comb$Row.names
comb<-comb[,-1]

vol_data<-data.frame(matrix(NA,(length(colnames(comb))-3),2))
rownames(vol_data)<-colnames(comb)[1:(length(colnames(comb))-3)]
colnames(vol_data)<-c("corr","pvalue")
for (i in 1:(length(colnames(comb))-3)){
  a<-cor.test(comb[,i],comb[,"mean_psi_clus2"],method = "pearson")
  vol_data[i,1] <-a$estimate
  vol_data[i,2] <-a$p.value
}
volcano<-subset(vol_data,select = c("corr","pvalue"))
threshold<-as.factor((abs(volcano$corr)>0.5) & volcano$pvalue<8.2e-5)   #其实最好是用bonferonni corrected P value
r03=ggplot(volcano,aes(corr,-log2(pvalue),colour=threshold))+geom_point()
r04=r03+labs(title="Volcanoplot of pearson corr with mean psi\n(Bonferroni-corrected test)")+theme(plot.title = element_text(hjust = 0.5))+xlim(-1,1)
r05=r04+geom_vline(xintercept=c(-0.5,0.5),linetype="dotted",size=1)+geom_hline(yintercept=-log2(8.2e-5),col="blue")

library(ggrepel)    #可以让字与字自动错开
label<-rownames(volcano)
label[which((abs(volcano$corr)<=0.4) | volcano$pvalue>=8.2e-5)]<-NA
r06=r05+geom_text_repel(aes(corr, -log2(pvalue), label = label))

########### 5.5 GO analysis ###########

library(topGO)   #画GO图用的

library(clusterProfiler)
library(Rgraphviz)
library(pathview)

#获得基因 symbol ID
#将symbolID转换成ENTREZID
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = nega_gene,
                       keytype = "SYMBOL",
                       column = "ENTREZID")   
#有67个基因没有转换成功
DEG.entrez_id = na.omit(DEG.entrez_id)
#biological process 富集分析
#gene_mat2<-read.csv (file = "/Users/chenxinyi/Desktop/gene_mat2.csv",sep=" ",header = TRUE,row.names = 1,stringsAsFactors = FALSE)
ref.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = colnames(var_expr),    #以那605个基因作为背景
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
                           qvalueCutoff = 0.5,
                           universe = ref.entrez_id)
## 画图
barplot(erich.go.CC_ref)

erich.go.MF_ref = enrichGO(gene = DEG.entrez_id,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID",
                           ont = "MF",
                           pvalueCutoff = 0.5,
                           qvalueCutoff = 0.5,
                           universe = ref.entrez_id)  # universe: background genes
barplot(erich.go.MF_ref)

erich.KEGG_ref=enrichKEGG(gene=DEG.entrez_id,
                          organism = "hsa",
                          keyType = "kegg",
                          pvalueCutoff = 0.5,
                          qvalueCutoff = 0.5,
                          universe = ref.entrez_id)
barplot(erich.KEGG_ref)


######### 5.6 FInd pathways correlated with global AS ##########
#先看iPSC
comb2<-merge(nes_Hs,all_psi_new,by=0,all=FALSE)
rownames(comb2)<-comb2$Row.names
comb2<-comb2[,-1]
comb2<-comb2[which(comb2$cell_type=="MN"),]

vol_data<-data.frame(matrix(NA,(length(colnames(comb2))-3),2))
rownames(vol_data)<-colnames(comb2)[1:(length(colnames(comb2))-3)]
colnames(vol_data)<-c("corr","pvalue")
for (i in 1:(length(colnames(comb2))-3)){
  a<-cor.test(comb2[,i],comb2[,"glbl_mean_psi"],method = "pearson")
  vol_data[i,1] <-a$estimate
  vol_data[i,2] <-a$p.value
}
volcano<-subset(vol_data,select = c("corr","pvalue"))
threshold<-as.factor((abs(volcano$corr)>0.5) & volcano$pvalue<1.7e-5)   #其实最好是用bonferonni corrected P value
r03=ggplot(volcano,aes(corr,-log2(pvalue),colour=threshold))+geom_point()
r04=r03+labs(title="Volcanoplot of pearson corr with mean psi in MN\n(Bonferroni-corrected test)")+theme(plot.title = element_text(hjust = 0.5))+xlim(-1,1)
r05=r04+geom_vline(xintercept=c(-0.5,0.5),linetype="dotted",size=1)+geom_hline(yintercept=-log2(1.7e-5),col="blue")

library(ggrepel)    #可以让字与字自动错开
label<-rownames(volcano)
label[which((abs(volcano$corr)<=0.5) | volcano$pvalue>=1.7e-5)]<-NA
r06=r05+geom_text_repel(aes(corr, -log2(pvalue), label = label))

posi_top<-rownames(volcano)[which(volcano$corr>0.45)] 
#720条pathway，都不属于Dpath,都不属于之前找的highly variable genes
length(setdiff(posi_top,Dpath))
posi_top %in% VariableFeatures(sc_activ)

nega_top<-rownames(volcano)[which(abs(volcano$corr)>0.45 & volcano$corr<0 )] 
#651条pathway，都不属于Dpath
length(setdiff(nega_top,Dpath))

posi_top<-posi_top[order(posi_top,decreasing=TRUE)]  #进行排序
nega_top<-nega_top[order(nega_top,decreasing=FALSE)]  #进行排序


layout(matrix(1:3, 1, 3,byrow=T))
plotfunc(comb2,"GO_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM","protein localization to ER")
plotfunc(comb2,"GO_POSITIVE_REGULATION_OF_CYTOKINESIS","positive regulation of cytokinesis")
plotfunc(comb2,"GO_ACETYL_COA_METABOLIC_PROCESS","acetyl CoA metabolic process")


########### 6. plot the number of reads in all three cells ############

colors <- c("green","orange","brown")
cellname <- c("MN","NPC","iPSC")
regions <- c("junc13","junc12","junc23")

# Create the matrix of the values.
Values <- matrix(c(314827.8, 228192.4, 859024.4, 678984.3, 520195.4, 2394778, 806892.7, 623268.2, 2758233),nrow=3,ncol=3,byrow=TRUE)
# Give the chart file a name.
#png(file = "barchart_stacked.png")
# Create the bar chart.
barplot(Values,names.arg=cellname,xlab="cell type",ylab="sum of junction reads",col=colors)
# Add the legend to the chart.
legend("topleft", regions, cex=1.3, fill=colors)
 # Save the file.
