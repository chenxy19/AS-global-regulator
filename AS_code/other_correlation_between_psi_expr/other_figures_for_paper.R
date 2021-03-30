#sce mouseSC
expr_psi<-read.csv("/Users/chenxinyi/Desktop/output/data3_mouseSC/local_mouseSC18.csv",sep=" ")
expr_psi2<-expr_psi[,-ncol(expr_psi)]
mouseSC_18W_GENIE3_rank <-run_GENIE3(t(expr_psi2),"glbl_mean_psi","/Users/chenxinyi/Desktop/output/data3_mouseSC/mouseSC_18W_GENIE3_rank.csv")

norm_expr<-read.csv("/Users/chenxinyi/Desktop/output/data4_mouseSC2/mouseSC2_sce_expr_qNSC1.csv",sep=" ")
mouseSC2_psi_bimodal<-read.csv('/Users/chenxinyi/Desktop/output/data4_mouseSC2/mouseSC2_psi_bimodal.csv',sep=" ",stringsAsFactors = F)
mouseSC2_expr<-as.data.frame(t(norm_expr))
mouseSC2_expr_psi <- merge_by_rownames(mouseSC2_expr,mouseSC2_psi_bimodal)
mouseSC2_expr_psi2<-mouseSC2_expr_psi[,-ncol(mouseSC2_expr_psi)]
mouseSC2_GENIE3_rank <-run_GENIE3(t(mouseSC2_expr_psi2),"glbl_mean_psi","/Users/chenxinyi/Desktop/output/data4_mouseSC2/mouseSC2_qNSC1_GENIE3_rank.csv")

corr_volcano<-find_gene_high_correlation(mouseSC2_expr_psi,"glbl_mean_psi")
corr_rank<-rownames(corr_volcano)[order(corr_volcano$corr,decreasing=TRUE)]  #进行排序
write.table (corr_rank, file = "/Users/chenxinyi/Desktop/output/data4_mouseSC2/mouseSC2_qNSC1_corr_rank.csv", row.names =TRUE, col.names =TRUE, quote =TRUE)


######### sup fig. distribution of read counts #########

data1<-read.csv("/Users/chenxinyi/Downloads/filereport_read_run_PRJNA339740_tsv-2.txt",sep="\t")
data2<-read.csv("/Users/chenxinyi/Downloads/filereport_read_run_PRJNA314949_tsv-3.txt",sep="\t")
data3<-read.csv("/Users/chenxinyi/Downloads/filereport_read_run_PRJNA345438_tsv-5.txt",sep="\t")
data4<-read.csv("/Users/chenxinyi/Downloads/filereport_read_run_PRJNA475579_tsv-6.txt",sep="\t")
data5<-read.csv("/Users/chenxinyi/Downloads/filereport_read_run_PRJEB26646_tsv-4.txt",sep="\t")
#可以draw violin plot

median(data1[grepl("^MN single",data1$sample_title,),"read_count"])
library(dplyr)
datasum=group_by(data2[,c("sample_accession","base_count","read_count")], sample_accession) %>% summarize_each(funs(sum))
median(datasum$read_count)

datasum4=group_by(data4[grep("single_cell",data4$sample_title),c("sample_accession","base_count","read_count")], sample_accession) %>% summarize_each(funs(sum))
median(datasum4$read_count)

datasum5=group_by(data5[,c("sample_accession","read_count","base_count")], sample_accession) %>% summarize_each(funs(sum))
median(datasum5$base_count)

data<-data1[,c("read_count","sample_title")]
data$sample_label<-0
data$sample_label[grepl("^MN single",data1$sample_title,)]<-"1_hMN"
data$sample_label[grepl("^NPC single",data1$sample_title,)]<-"1_hNPC"
data$sample_label[grepl("^iPSC single",data1$sample_title,)]<-"1_hiPSC"
#data<-data[,c("read_count","sample_label")]

data_add<-datasum
data_add$sample_label<-c("2_hHSC")
data_add<-data_add[,c("read_count","sample_accession","sample_label")]
colnames(data_add)<-c("read_count","sample_title","sample_label")
data<-rbind(data,data_add)

data_add<-data3
data_add$sample_label<-c("3_mHSC")
data_add<-data_add[,c("read_count","sample_accession","sample_label")]
colnames(data_add)<-c("read_count","sample_title","sample_label")
data<-rbind(data,data_add)

data_add<-datasum4
data_add$sample_label<-c("4_mNSC")
data_add<-data_add[,c("read_count","sample_accession","sample_label")]
colnames(data_add)<-c("read_count","sample_title","sample_label")
data<-rbind(data,data_add)

data_add<-datasum5
data_add$sample_label<-c("5_hESC")
data_add<-data_add[,c("read_count","sample_accession","sample_label")]
colnames(data_add)<-c("read_count","sample_title","sample_label")
data<-rbind(data,data_add)

library(ggplot2) 
ggplot(data[which(data$sample_label!=0),], aes(x = sample_label , y = read_count , fill = sample_label )) +
  geom_violin(alpha = 0.5,scale="width",aes(linetype=NA)) + 
  geom_jitter(shape=21,aes(fill=sample_label),position = position_jitter(width = 0.2))+
  xlab("Datasets")+ylab("Average read counts per cell")+
  theme_bw()+theme(legend.position = "none")

 
######### sup fig. dcorrelation between mean psi #########
#统一画吧
#file="/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/iPSC_2p_2meth.csv"
files=c("/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/MN_2p_2meth.csv",
        "/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/NPC_2p_2meth.csv",
        "/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/iPSC_2p_2meth.csv",
        "/Users/chenxinyi/Desktop/output/data2_hHSC/HSC_2p_2meth.csv",
        "/Users/chenxinyi/Desktop/output/data5_human3/human3_2p_2meth.csv",
        "/Users/chenxinyi/Desktop/output/data3_mouseSC/mouseSC_2p_2meth.csv",
        "/Users/chenxinyi/Desktop/output/data4_mouseSC2/mouseSC2_2p_2meth.csv")
cell_type=c("1_hMN","1_hNPC","1_hiPSC","2_hHSC","5_hESC","3_mHSC","4_mNSC")
layout(matrix(1:1, 1, 1,byrow=T))
for (i in 1:7){
  iPSC2p_2meth <- read.csv(files[i],header = TRUE,stringsAsFactors = FALSE, sep=" ") 
  head(iPSC2p_2meth)
  iPSC2p_2meth<-iPSC2p_2meth[which(!is.na(iPSC2p_2meth$bi_global.aver)),]
  colnames(iPSC2p_2meth )<-c("bi_global-aver","bi_gene-aver","all_global-aver","all_gene-aver")
  library(corrplot)
  splicr2_cor<-cor(iPSC2p_2meth,method = "pearson") #计算Pearson相关系数矩阵（用于连续性变量，变量服从正态分布）
  col3 <- colorRampPalette(c("blue", "white", "red")) 
  
  res1 <- cor.mtest(iPSC2p_2meth, conf.level = .95)
  corrplot(splicr2_cor, p.mat = res1$p, col=col3(10),insig = "label_sig",type="upper",tl.pos="left",tl.cex = 0.75,
           sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white",title=cell_type[i],mar = c(1,1,1,1))
  cor.plot <- corrplot(corr = splicr2_cor,add=TRUE, type="lower",col=col3(10),method="color",addCoef.col="black",diag=FALSE,tl.pos="n", cl.pos="n",number.cex = 0.7)
  
}


#不同数据集的mean psi的分布也值得一画
#看一下各个数据集之间bimodal AS 有多少重叠（这个可以本地做）
#每个细胞里算的bimodal gene都不一样，这是主要问题！！！应该看一下每个单细胞里能检测到的bimodal events是不是比较固定的
#如果是用gene average的话，最好是固定一群基因
# plot how sparse the psi data are
#liana的paper可能真的是很重要的
#只把NMD和其他几种的相关性求一下，其他都忽略


dir="/Users/chenxinyi/Desktop/output/"
prefix="data1_MN+NPC+iPSC/MN_"
bimodal_file = paste(dir,prefix,"se_bi_psi.csv",sep="")
bimodal= read.csv(bimodal_file,header = TRUE,stringsAsFactors = FALSE) 

event_file=paste(prefix,"events.csv",sep="")
events= read.csv(event_file,header = TRUE,stringsAsFactors = FALSE) 

######### sup fig. distribution of mean psi #########
all_psi=read.csv(files[1],header = TRUE,stringsAsFactors = FALSE, sep=" ") 
all_psi<-all_psi[which(!is.na(all_psi[,1])),]
colnames(all_psi )<-c("bi_global_aver","bi_gene_aver","all_global_aver","all_gene_aver")
all_psi$label=cell_type[1]

for (i in 2:7){
  iPSC2p_2meth <- read.csv(files[i],header = TRUE,stringsAsFactors = FALSE, sep=" ") 
  iPSC2p_2meth<-iPSC2p_2meth[which(!is.na(iPSC2p_2meth[,1])),]
  colnames(iPSC2p_2meth)<-c("bi_global_aver","bi_gene_aver","all_global_aver","all_gene_aver")
  iPSC2p_2meth$label=cell_type[i]
  all_psi <- rbind(all_psi , iPSC2p_2meth)
}

library(reshape2)
all_psi_melt<-melt(all_psi)
all_psi_melt$whole_label<-all_psi_melt$label+all_psi_melt$variable

library(ggplot2) 
#分组的violinplot，但是没有了
ggplot(all_psi_melt, aes(x = label , y = value , fill = label )) +
  geom_violin(alpha = 0.5,scale="width",aes(fill = variable), trim = FALSE) + 
  #geom_jitter(shape=21,aes(fill=label),position = position_jitter(width = 0.2))+
  xlab("Datasets")+ylab("Average global-mean psi for bimodal genes")+
  theme_bw()+ labs(fill = "Type of mean psi")#+theme(legend.position = "none")
   

############# permutation of psi file and recalculate mean psi ###########
# random permutation for each AS events by shuffle row-wise

MN_psi= read.csv("/Users/chenxinyi/Desktop/output/data5_human3/psi.csv",sep=",")
rownames(MN_psi)<-MN_psi$sample_id
MN_psi<-MN_psi[,-1]
mean_psi_all<-apply(MN_psi, 1, mean_non_na)
for (i in 1:200){
  MN_psi_rand = lapply(MN_psi, function(x) { sample(x) })  #每一列分别去打乱，没有coupling
  MN_psi_rand<-as.data.frame(MN_psi_rand)
  mean_psi_perm<-apply(MN_psi_rand, 1, mean_non_na)
  mean_psi_all<-cbind(mean_psi_all, mean_psi_perm)
}

mean_psi_all<-as.data.frame(mean_psi_all)
colnames(mean_psi_all)[2:ncol(mean_psi_all)]<-paste("mean_psi",1:200,sep="_")

#去掉outlier再作图
obs<-mean_psi_all$mean_psi_all[!is.na(mean_psi_all$mean_psi_all)]
lower_bound <- median(obs) - 3 * mad(obs)
upper_bound <- median(obs) + 3 * mad(obs)
obs_test<-obs[which(obs<upper_bound & obs>lower_bound)]

hist(melt(mean_psi_all[,2:200])$value,freq=FALSE,breaks = 10,xlim=range(0.6,0.9),col=rgb(1,0,0,1/4),main = "Observed & permutated mean psi (5_hESC)",xlab = "mean_psi")
hist(obs_test,freq=FALSE,breaks = 10,add=T,col=rgb(0,0,1,1/4),xlab = "mean_psi",xlim=range(0.56,0.72))
legend("topright", c("Permutated (200)", "Observed"),fill=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)),bty="n")  #bty="n"是legend周围的边框不要

# KS test
ks.test(melt(mean_psi_all[,2:200])$value, obs_test )
# shapiro.test
shapiro.test(obs_test)  # 一般都非正态
# F test（需要假设正态）
#var.test(melt(mean_psi_all[,2:17])$value, mean_psi_all$mean_psi_all, alternative = "two.sided")
# Levene test
library(car)

all_melt<-melt(mean_psi_all[,2:200])
all_melt$variable<-"mean_psi_perm"
obs_test<-melt(obs_test)
obs_test$variable<-"mean_psi_obs"
all_melt<-all_melt[,c("value","variable")]
all_melt<-rbind(obs_test,all_melt)

leveneTest(all_melt$value, all_melt$variable)
#前者是值，后者是分组，需要同长

#Two-sample Kolmogorov-Smirnov test
#data:  melt(mean_psi_all[, 2:20])$value and mean_psi_all$mean_psi_all
#D = 0.24763, p-value = 0.001639  (是显著不同的)
#alternative hypothesis: two-sided

############## NMD mean psi ###############
MN_psi= read.csv("/Users/chenxinyi/Desktop/output/NMD_psi.csv",sep=",")
rownames(MN_psi)<-MN_psi$sample_id
MN_psi<-MN_psi[,-1]
mean_psi_all<-apply(MN_psi, 1, mean_non_na)


