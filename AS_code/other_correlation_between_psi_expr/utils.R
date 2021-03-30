"""
Useful functions.
"""

library(reshape2)
library(ggplot2)
library(Seurat)
library(ggpmisc)
library(philentropy)

#random forest identify the top influential variables
library(GENIE3)

#GO
require(org.Hs.eg.db)  #这个包里存有人的注释文件
library(topGO)   #画GO图用的
library(clusterProfiler)
library(Rgraphviz)
library(pathview)

#Heatmap
library(gplots)
library(ComplexHeatmap)

# gene set analysis
# usage in single_cell_pathway_activity.R
library(viper)    
library(msigdbr)

calculate_psi<-function(prefix,mat_junc_cell,SRRlist=NULL){    #函数可以带有缺省值或者默认值
  junc_reads_file=paste("/Users/chenxinyi/Desktop/output/",prefix,"_mean_psi/reads.csv",sep = "")
  junc_reads= read.csv(junc_reads_file,header = TRUE,stringsAsFactors = FALSE) 
  if (is.null(SRRlist)){
    SRRlist<-unique(junc_reads$sample_id)
  }
  
  event_file=paste("/Users/chenxinyi/Desktop/output/",prefix,"_mean_psi/events.csv",sep = "")
  events= read.csv(event_file,header = TRUE,stringsAsFactors = FALSE) 
  junctions<-events[,c("event_id","junction13","junction12","junction23")]    #其实unique只有8690,which is reasonable
  junctions<-unique(junctions)
  junctions$event_id<- gsub("[^0-9A-z]",".",junctions$event_id)
  ASlist<-junctions$event_id
  
  psi_junc_cell<-data.frame(matrix(NA,length(ASlist),length(SRRlist)))
  colnames(psi_junc_cell)<-SRRlist
  rownames(psi_junc_cell)<-ASlist
  for (j in 1:length(SRRlist)){
    for (i in 1:length(ASlist)){
      junc13<-as.numeric(mat_junc_cell[3*i-2,j])
      junc12<-as.numeric(mat_junc_cell[3*i-1,j])
      junc23<-as.numeric(mat_junc_cell[3*i,j])
      if (junc13 + junc12 + junc23  >= 10){
        psi_junc_cell[i,j]<-(junc12 + junc23)/(junc12 + junc23 + 2*junc13)
      }
    }
  }
  return(psi_junc_cell)
}

calculate_global_mean_psi<-function(prefix,mat_junc_cell,SRRlist=NULL){    #函数可以带有缺省值或者默认值
  junc_reads_file=paste("/Users/chenxinyi/Desktop/output/",prefix,"_mean_psi/reads.csv",sep = "")
  junc_reads= read.csv(junc_reads_file,header = TRUE,stringsAsFactors = FALSE) 
  if (is.null(SRRlist)){
    SRRlist<-unique(junc_reads$sample_id)
  }
  
  event_file=paste("/Users/chenxinyi/Desktop/output/",prefix,"_mean_psi/events.csv",sep = "")
  events= read.csv(event_file,header = TRUE,stringsAsFactors = FALSE) 
  junctions<-events[,c("event_id","junction13","junction12","junction23")]    #其实unique只有8690,which is reasonable
  junctions<-unique(junctions)
  junctions$event_id<- gsub("[^0-9A-z]",".",junctions$event_id)
  ASlist<-junctions$event_id
  
  sum_junc<-data.frame(matrix(0,3,length(SRRlist)))
  colnames(sum_junc)<-SRRlist
  rownames(sum_junc)<-c("junc13","junc12","junc23")
  for (j in 1:length(SRRlist)){
    for (i in 1:length(ASlist)){
      junc13<-as.numeric(mat_junc_cell[3*i-2,j])
      junc12<-as.numeric(mat_junc_cell[3*i-1,j])
      junc23<-as.numeric(mat_junc_cell[3*i,j])
      if (junc13 + junc12 + junc23  >= 10){
        sum_junc[1,j] <- sum_junc[1,j] + junc13
        sum_junc[2,j] <- sum_junc[2,j] + junc12
        sum_junc[3,j] <- sum_junc[3,j] + junc23
      }
    }
  }
  glbl_mean_psi<-apply(sum_junc,2,divide)
  
  return(glbl_mean_psi)
}

divide<-function(a){
  psi<-(a[2]+a[3])/(a[2]+a[3]+2*a[1])
  return(psi)
}


downsample<-function(MN_junc_mat){
  MN_junc_mat_QC <- MN_junc_mat[,which(colnames(MN_junc_mat) %in% filtered$x)]
  prop1<-rep(NA,times=ncol(MN_junc_mat_QC))
  for (i in 1:length(prop1)){
    prop1[i] <- 165703/sum(MN_junc_mat_QC[,i])
  }
  MN_junc_mat_QC_downsampled<-downsampleMatrix(as.matrix(MN_junc_mat_QC),prop=prop1 ,bycol=TRUE)   #downSample的函数只接受matrix，不接受dataframe!
  return(MN_junc_mat_QC_downsampled)
}

merge_by_rownames<-function(dt1,dt2,all_choice=FALSE){
  a<-merge(dt1,dt2,by=0,all=all_choice)
  rownames(a)<-a[,1]
  a<-a[,-1]
  return(a)
}

sum_non_na<-function(a){
  num_non_na<-sum(!is.na(a))
  return(num_non_na)
}

mean_non_na<-function(a){
  return(mean(a, na.rm=TRUE))
}

#可以用apply(X,1,mean_na.rm)代替
calculate_gene_average_psi<-function(file){
  se_psi= read.csv(file,header = TRUE,row.names = 1,stringsAsFactors = FALSE) 
  #se_psi<-se_psi[!rownames(se_psi) %in% c("SRR4047260", "SRR4047278","SRR4047295","SRR4047315","SRR4047390","SRR4047332","SRR4047331", "SRR4047323"),] 
  se_psi_mean<-se_psi
  se_psi_mean$mean_psi<-0
  for (i in 1:nrow(se_psi)){
    se_psi_mean$mean_psi[i]<-mean(as.numeric(se_psi[i,which(!is.na(se_psi[i,]))]))   #na.rm=TRUE
  }
  title=paste("Histogram of mean psi across single ",prefix,sep="")
  hist(se_psi_mean$mean_psi,breaks = 20,xlab = "mean_psi",main = title)
  return(se_psi_mean)
}

merge_into_all_matrix<-function(prefix,se_psi_mean){
  se_psi_mean$cell_type<-prefix
  a<-merge(all_globl_mean_psi,se_psi_mean,by=0,all=FALSE)
  rownames(a)<-a$Row.names
  a<-a[,c("glbl_mean_psi","mean_psi","cell_type" )]
  colnames(a)[2]<-"mean_psi_genes"
  return(a)
}

find_gene_high_correlation<-function(comb,colname_mean_psi,
                                     cor_method="pearson",
                                     cor_threshold=0,
                                     tail_col=2){
  vol_data<-data.frame(matrix(NA,(ncol(comb)-tail_col),2))
  rownames(vol_data)<-colnames(comb)[1:(ncol(comb)-tail_col)]
  colnames(vol_data)<-c("corr","pvalue")
  for (i in 1:(ncol(comb)-tail_col)){
    a<-cor.test(comb[,i],comb[,colname_mean_psi],method = cor_method)
    vol_data[i,1] <-a$estimate
    vol_data[i,2] <-a$p.value
  }
  volcano<-subset(vol_data,select = c("corr","pvalue"))
  n_Bonferroni_correction=ncol(comb)-tail_col
  pvalue_threshold<-0.05/n_Bonferroni_correction
  threshold<-as.factor((abs(volcano$corr)>cor_threshold) & volcano$pvalue<pvalue_threshold)   #其实最好是用bonferonni corrected P value
  r03=ggplot(volcano,aes(corr,-log2(pvalue),colour=threshold))+geom_point()
  r04=r03+labs(title=paste("Volcanoplot of ",cor_method," corr with mean psi\n(Bonferroni-corrected test)",sep=""))+theme(plot.title = element_text(hjust = cor_threshold))+xlim(-1,1)
  r05=r04+geom_vline(xintercept=c(-1*cor_threshold,cor_threshold),linetype="dotted",size=1)+geom_hline(yintercept=-log2(pvalue_threshold),col="blue")
  
  return(volcano)
  
  library(ggrepel)    #可以让字与字自动错开
  label<-rownames(volcano)
  label[which((abs(volcano$corr)<=cor_threshold) | volcano$pvalue>=pvalue_threshold)]<-NA
  r06=r05+geom_text_repel(aes(corr, -log2(pvalue), label = label))
  r06
  
  #posi_gene<-label[which(volcano$corr>cor_threshold & volcano$pvalue<pvalue_threshold)]
  #nega_gene<-label[which(abs(volcano$corr)>cor_threshold & volcano$corr<0 & volcano$pvalue<pvalue_threshold)]
  #e <- list(posi_gene,nega_gene,r06)  #因为函数只能返回一个值，所以给打包成一个对象
  #return(e)
}

find_gene_corr_FDR<-function(comb,colname_mean_psi,
                                     cor_method="pearson",
                                     cor_threshold=0,
                             tail_col=3){
  p_threshold=0.05
  vol_data<-data.frame(matrix(NA,(ncol(comb)-tail_col),2))
  rownames(vol_data)<-colnames(comb)[1:(ncol(comb)-tail_col)]
  colnames(vol_data)<-c("corr","pvalue")
  for (i in 1:(ncol(comb)-tail_col)){
    a<-cor.test(comb[,i],comb[,colname_mean_psi],method = cor_method)
    vol_data[i,1] <-a$estimate
    vol_data[i,2] <-a$p.value
  }
  volcano<-subset(vol_data,select = c("corr","pvalue"))
  volcano$pvalue_BH <- p.adjust(volcano$pvalue, "BH")
  threshold<-as.factor(volcano$pvalue_BH < p_threshold)   
  r03=ggplot(volcano,aes(corr,-log2(pvalue_BH),colour=threshold))+geom_point()
  r04=r03+labs(title=paste("Volcanoplot of ",cor_method," corr with mean psi\n(BH corrected)",sep=""))+theme(plot.title = element_text(hjust = p_threshold))+xlim(-1,1)
  #r05=r04+geom_vline(xintercept=c(-1*cor_threshold,cor_threshold),linetype="dotted",size=1)+geom_hline(yintercept=-log2(pvalue_threshold),col="blue")
  
  library(ggrepel)    #可以让字与字自动错开
  label<-rownames(volcano)
  label[which(volcano$pvalue_BH >= p_threshold)]<-NA
  r05=r04+geom_text_repel(aes(corr, -log2(pvalue_BH), label = label))
  
  posi_gene<-label[which(volcano$corr>cor_threshold & volcano$pvalue_BH < p_threshold)]
  nega_gene<-label[which(abs(volcano$corr)>cor_threshold & volcano$corr<0 & volcano$pvalue_BH < p_threshold)]
  e <- list(posi_gene,nega_gene,r05)  #因为函数只能返回一个值，所以给打包成一个对象
  return(e)
}

plotfunc<-function(comb,GO,x){
  a<-cor.test(comb[,GO],comb$glbl_mean_psi)
  print(a)
  plot(comb[,GO],comb$glbl_mean_psi,xlab=x,ylab="mean_psi (by sum of junctions)",main=paste("Pearson corr=",round(a$estimate,5), "\np-value=",round(a$p.value,5),sep=""))
  LM<-lm(comb$glbl_mean_psi ~ comb[,GO])
  abline(LM)
}

calculate_dist_with_na<-function(a,b,dist_method="euclidean"){
  nonna<-which(!is.na(a) & !is.na(b))
  vect<-data.frame(matrix(NA,2,length(nonna)))
  vect[1,]<-a[nonna]
  vect[2,]<-b[nonna]
  vect_dist<-dist(vect,method=dist_method)**2/length(nonna)    #平均euclidean distance
  return(vect_dist)
}

is_nega<-function(a){
  is_neg=""
  if (a<0){
    is_neg="-"
  }
  return(is_neg)
}

gene_mean_psi_scatter<-function(gene){
  row_number<-which(colnames(all_SF_psi)==gene)
  cor_data<-data.frame(NA,2,3)
  i=1
  for (cell in c("MN","NPC","iPSC")){
    selected_lines<-which(all_SF_psi$cell_type==cell)
    cor_test<-cor.test(all_SF_psi[selected_lines,row_number] , all_SF_psi$glbl_mean_psi[selected_lines], method = "pearson")
    cor_data[1,i]<-cor_test$estimate
    cor_data[2,i]<-cor_test$p.value
    i=i+1
  }
  
  labels=paste("MN_cor=",is_nega(cor_data[1,1]),abs(round(cor_data[1,1],5)) ,",p=",round(cor_data[2,1],5),
               "\nNPC_cor=",is_nega(cor_data[1,2]),abs(round(cor_data[1,2],5)) ,",p=",round(cor_data[2,2],5),
               "\niPSC_cor=",is_nega(cor_data[1,3]),abs(round(cor_data[1,3],5)) ,",p=",round(cor_data[2,3],5),
               sep="")
  
  pic<-ggplot(all_SF_psi,aes(x=all_SF_psi[,row_number],y=glbl_mean_psi,color=as.factor(cell_type))) + geom_point(shape=19) +
    geom_point() + 
    scale_colour_brewer(palette = "Set1") +
    geom_smooth(method="lm",show.legend = TRUE,level=0.95) +
    labs(title=labels, x=gene, y="mean psi (by summing up junction reads)")
  #scale_fill_discrete(breaks=c("5_n1","5_n2","3_n1","3_n2")) +
  return(pic)
}

plot_venn<-function(listA,listB,listC,filename="plot_venn.tiff"){
  library(VennDiagram)
  venn.plot <- venn.diagram(
    x = list(
      differential_genes = listA ,
      negative_cor_psi1 = listB ,
      negative_cor_psi2= listC
    ),
    filename = filename,
    col = "transparent",
    fill = c("red", "blue", "green"),
    alpha = 0.5,
    label.col = c("darkred", "white", "darkblue", "white",
                  "white", "white", "darkgreen"),
    cex = 2.5,#内标签的字体大小
    fontfamily = "serif",
    fontface = "bold",
    cat.default.pos = "outer",#设置标签在圆外面
    #cat.col = c("darkred", "darkblue", "darkgreen"),
    cat.cex = 1.5,#外标签的字体大小
    cat.fontfamily = "serif",
    cat.dist = c(0.05, 0.05, 0.05),#相对圆圈的位置
    cat.pos = c(-20,20,180)  #相对12点方向旋转的角度
  )
}

GO_analysis<-function(listA, background){
  require(org.Hs.eg.db)  
  library(topGO)   
  library(clusterProfiler)
  library(Rgraphviz)
  library(pathview)
  DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                         keys = listA,
                         keytype = "SYMBOL",
                         column = "ENTREZID")   
  no_NA <- sum(is.na(DEG.entrez_id))
  DEG.entrez_id = na.omit(DEG.entrez_id)

  ref.entrez_id = mapIds(x = org.Hs.eg.db,
                         keys = background,
                         keytype = "SYMBOL",
                         column = "ENTREZID")  
  no_NA_bg <- sum(is.na(ref.entrez_id))
  ref.entrez_id = na.omit(ref.entrez_id)
  
  erich.go.BP_ref = enrichGO(gene = DEG.entrez_id,
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID",
                             ont = "BP",
                             pvalueCutoff = 0.5,
                             qvalueCutoff = 0.5,
                             universe = ref.entrez_id)  # universe: background genes
  plotBP<- dotplot(erich.go.BP_ref)
  
  erich.go.CC_ref = enrichGO(gene = DEG.entrez_id,
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID",
                             ont = "CC",
                             pvalueCutoff = 0.5,
                             qvalueCutoff = 0.5,
                             universe = ref.entrez_id)
  plotCC <- barplot(erich.go.CC_ref)
  
  erich.go.MF_ref = enrichGO(gene = DEG.entrez_id,
                             OrgDb = org.Hs.eg.db,
                             keyType = "ENTREZID",
                             ont = "MF",
                             pvalueCutoff = 0.5,
                             qvalueCutoff = 0.5,
                             universe = ref.entrez_id)  # universe: background genes
  plotMF<- barplot(erich.go.MF_ref)
  
  erich.KEGG_ref=enrichKEGG(gene=DEG.entrez_id,
                            organism = "hsa",
                            keyType = "kegg",
                            pvalueCutoff = 0.5,
                            qvalueCutoff = 0.5,
                            universe = ref.entrez_id)
  plotKEGG<- barplot(erich.KEGG_ref)
  return(list(plotBP, plotCC, plotMF, plotKEGG, no_NA ,no_NA_bg ))
}

run_GENIE3 <- function(expr_psi_mat, target_name,file_name) {   #matrix should be gene*sample
  expr_psi_mat[is.na(expr_psi_mat)] <- 0
  weightMat<-GENIE3(as.matrix(expr_psi_mat),targets=target_name,nCores=1,nTrees=10000)  #1000是default
  weight_index<-sort(weightMat, decreasing = TRUE, index.return=TRUE)
  print(rownames(expr_psi_mat)[weight_index$ix[1:50]])
  top_gene<-rownames(expr_psi_mat)[weight_index$ix]
  write.table (top_gene, file = file_name, row.names =TRUE, col.names =TRUE, quote =TRUE)
  return(top_gene)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
