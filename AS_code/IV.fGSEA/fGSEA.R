"""
Run fast GSEA to find gene sets wth high correlation with mean psi. Fast GSEA take uneven steps according to the enrichment score.

   Input: Zscore.csv

"""


library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)

m_df<- msigdbr(species = "Homo sapiens", category = "C5")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

prefixes=c("MN","NPC","iPSC","hHSC","mHSC","mNSC","hESC")
for (prefix in prefixes){
  zscore_corr <- read.csv(paste("/Users/chenxinyi/Desktop/output/permt_zscore/",prefix,"_Zscore.csv",sep=""))
  rownames(zscore_corr) <- zscore_corr$X
  
  path_ES = c()
  
  for (n_col in 2:ncol(zscore_corr)){

    genes_rank <- zscore_corr[,c(1,n_col)]
    colnames(genes_rank)<-c("gene","zscore")
    
    genes_rank<-genes_rank %>%
      arrange(desc(zscore)) 
      
    ranks<- tibble::deframe(genes_rank)
    
    # gene sets too big or to small is meaningless
    #eps=0 allows accurate estimation.
    # when nperm (permutation) is not called, multilevel is run.
    fgseaRes<- fgsea(fgsea_sets, stats = ranks, minSize=15, maxSize=200,eps=0)
    
    fgseaResTidy <- fgseaRes %>%
      as_data_frame() %>%
      arrange(desc(NES))
    
    # NES为每个基因子集ES根据基因集大小标准化的富集分数，ES反应基因集在排序列表两端的富集程度
    # p.adjust：'BH' 校准后的P值
    
    if (is.null(nrow(path_ES))){
      path_ES<-data.frame(fgseaResTidy[,c("pathway", "NES")])
      colnames(path_ES)<-c("pathway", "observed")
    }else{
      path_ES<-merge(path_ES, data.frame(fgseaResTidy[,c("pathway", "NES")]), by="pathway", all=TRUE)
      colnames(path_ES)[ncol(path_ES)]<-paste("iter", n_col-2,sep="")
    }
    write.table(path_ES, file = paste("/Users/chenxinyi/Desktop/output/permt_zscore/", prefix, "_pathwayNES.csv",sep=""))
  }
}


obs_rank<-function(list){
  return(rank(list)["observed"])
}

for (prefix in prefixes){
  NES <- read.csv(paste("/Users/chenxinyi/Desktop/output/permt_zscore/", prefix, "_pathwayNES.csv",sep=""), sep=" ",  row.names="pathway")
  NES<-NES[,-1]
  
  obs_rank_path<- apply(NES, 1, obs_rank)
  obs_rank_path_df <- as.data.frame(obs_rank_path)
  rank_500 <- rownames(obs_rank_path_df)[which(obs_rank_path_df$obs_rank_path>=495)]
  rank_1 <- rownames(obs_rank_path_df)[which(obs_rank_path_df$obs_rank_path<=5)]
  # 写进文件里
  write.table(rank_500, file = paste("/Users/chenxinyi/Desktop/output/permt_zscore/", prefix, "_rank500.csv",sep=""))
  write.table(rank_1, file = paste("/Users/chenxinyi/Desktop/output/permt_zscore/", prefix, "_rank1.csv",sep=""))
  
}


prefixes=c("MN","NPC","iPSC","hHSC","mouseSC","mouseSC2","human3")
for (prefix in prefixes){
  zscore_corr <- read.csv(paste("/Users/chenxinyi/Desktop/output/zscore/",prefix, "_Zscore.csv", sep=""))
  colnames(zscore_corr)<-c("feature","corr","zscore")
  
  # can use either corr or zscore
  genes_rank <- zscore_corr %>%
    arrange(desc(zscore)) %>%
    dplyr::select(feature, zscore)
  
  ranks<- tibble::deframe(genes_rank)
  
  # gene sets too big or to small is meaningless
  #eps=0 allows accurate estimation.
  # when nperm (permutation) is not called, multilevel is run.
  fgseaRes<- fgsea(fgsea_sets, stats = ranks, minSize=15, maxSize=200,eps=0)
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  #fgseaResTidy %>%
  #  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  #  arrange(padj) %>%
  #  head(n=20)
  
  pathways <- fgseaResTidy %>%
    dplyr::select(-leadingEdge, -ES) %>%
    arrange(padj)
  #write.table(pathways, file = paste("/Users/chenxinyi/Desktop/output/zscore/", prefix, "_pathway_gene.csv",sep=""))
  
  # need enough iteration to achieve significant result
  ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill= NES < 2)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste("Hallmark pathways NES from GSEA (",prefix,")",sep="")) +
    theme_minimal() 
  
  plot_a<-list()
  fontsize <- theme(axis.title=element_text(size=8), title=element_text(size=8))
  #title<-c("cytosolic ribosome", "cytosolic large ribosome\n subunit", "cotranslational protein target \nto membrane ",
  #         "protein localization to ER", "establishment of protein\n localization","mRNA NMD")
  title2<-c("polysomal ribosome", "protein localization\n to ER", "protein localization to ER",
            "mRNA NMD", "cotranslational protein target \nto membrane", "cytosolic ribosome")
  for (i in 1:6){
    plot_a[[i]]<-plotEnrichment(fgsea_sets[[pathways$pathway[i]]],
                                ranks) + labs(title=paste(pathways$pathway[i], " (", prefix," +)",sep="")) + fontsize
  }
  # run useful_function first
  multiplot(plot_a[[1]],plot_a[[2]],plot_a[[3]],plot_a[[4]],plot_a[[5]],plot_a[[6]], cols=2)
  
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  pdf(paste("/Users/chenxinyi/Desktop/output/zscore/", prefix, "_glbl.pdf",sep=""))
  plotGseaTable(fgsea_sets[topPathways], ranks, fgseaRes, 
                gseaParam = 0.5)
  dev.off()
}

# can get rid of repeating pathways 
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                      fgsea_sets, ranks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
dev.off()
plotGseaTable(fgsea_sets[mainPathways], ranks, fgseaRes, 
              gseaParam = 0.5,colwidths = c(10, 3, 0.8, 1.2, 1.2))

lead_gene<-fgseaRes[order(pval)][pathway %in% mainPathways]$leadingEdge
all_lead_gene<-c()
for (i in 1:length(lead_gene)){
  all_lead_gene<-union(all_lead_gene,lead_gene[[i]])
}


# only look at those that are < 0.05
# 66 genes for NPC, 3 collapsed pathway
# 416 genes for MN, 22 collapsed pathway
# 56 genes, 4 collapsed pathways


MN_all_lead_gene<-all_lead_gene
NPC_all_lead_gene<-all_lead_gene
mSC_all_lead_gene<-all_lead_gene
human3_all_lead_gene<-all_lead_gene
# Only 30 RPL and RPS are overlapping
conflict_gene<-intersect(MN_all_lead_gene, NPC_all_lead_gene)
NPC_zscore_corr<-zscore_corr
merged_zscore<-merge(MN_zscore_corr,NPC_zscore_corr,by="feature",suffixes = c(".MN",".NPC"), all=FALSE)
rownames(merged_zscore)<-merged_zscore$feature
plot(merged_zscore[conflict_gene,"zscore.MN"], merged_zscore[conflict_gene,"zscore.NPC"], main = "Overlapping leading genes", xlab="z score in MN",ylab="z score in NPC")
sele_zscore<-merged_zscore[union(MN_all_lead_gene, NPC_all_lead_gene),]
plot(sele_zscore[,"zscore.MN"], sele_zscore[,"zscore.NPC"], main = "Union of leading genes", xlab="z score in MN",ylab="z score in NPC")

sele_zscore$class<-0
sele_zscore[MN_all_lead_gene,"class"]<-"MN leading gene"
sele_zscore[NPC_all_lead_gene,"class"]<-"NPC leading gene"
sele_zscore[conflict_gene,"class"]<-"both"

ggplot(sele_zscore,aes(x=`zscore.MN`,y=`zscore.NPC`)) + geom_point(shape=19,alpha=0.6,aes(colour=as.factor(class))) +
  #scale_colour_brewer(palette = "Set1")+
  #scale_colour_gradientn(colours = terrain.colors(10)) +
  # geom_smooth(method="lm",show.legend = TRUE,level=0.95) +
  labs(title="Z score for leading genes")+
  xlab("Z score in MN") + ylab("Z score in NPC") +
  theme(plot.title = element_text(hjust = 0.5))
  

selected_GO_gene<-fgsea_sets$GO_CYTOSOLIC_RIBOSOME
sele_GO_zscore<-merged_zscore[selected_GO_gene,]
sele_GO_zscore$class<-"not leading gene"
sele_GO_zscore[MN_all_lead_gene,"class"]<-"MN leading gene"
sele_GO_zscore[NPC_all_lead_gene,"class"]<-"NPC leading gene"
sele_GO_zscore[conflict_gene,"class"]<-"both"
ggplot(sele_GO_zscore,aes(x=`zscore.MN`,y=`zscore.NPC`)) + geom_point(shape=19,alpha=0.8,aes(colour=as.factor(class))) +
  #scale_colour_brewer(palette = "Set1")+
  #scale_colour_gradientn(colours = terrain.colors(10)) +
  # geom_smooth(method="lm",show.legend = TRUE,level=0.95) +
  labs(title="Z score for leading genes GO_cytosolic_ribosome")+
  xlab("Z score in MN") + ylab("Z score in NPC") +
  theme(plot.title = element_text(hjust = 0.5))





conflict_gene<-intersect(MN_all_lead_gene, mSC_all_lead_gene)
merged_zscore<-merge(merged_zscore,zscore_corr,by="feature",suffixes = c("",".mSC"), all=FALSE)
colnames(merged_zscore)[6:7]<-c("corr.mSC", "zscore.mSC")
rownames(merged_zscore)<-merged_zscore$feature
sele_zscore<-merged_zscore[union(MN_all_lead_gene, mSC_all_lead_gene),]

sele_zscore$class<-0
sele_zscore[MN_all_lead_gene,"class"]<-"MN leading gene"
sele_zscore[mSC_all_lead_gene,"class"]<-"mSC leading gene"
sele_zscore[intersect(MN_all_lead_gene, mSC_all_lead_gene),"class"]<-"both"

ggplot(sele_zscore,aes(x=`zscore.MN`,y=`zscore.mSC`)) + geom_point(shape=19,alpha=0.6,aes(colour=as.factor(class))) +
  #scale_colour_brewer(palette = "Set1")+
  #scale_colour_gradientn(colours = terrain.colors(10)) +
  # geom_smooth(method="lm",show.legend = TRUE,level=0.95) +
  labs(title="Z score for leading genes")+
  xlab("Z score in MN") + ylab("Z score in mSC") +
  theme(plot.title = element_text(hjust = 0.5))

selected_GO_gene<-fgsea_sets$GO_CYTOSOLIC_RIBOSOME
sele_GO_zscore<-merged_zscore[selected_GO_gene,]
sele_GO_zscore$class<-"not leading gene"
sele_GO_zscore[MN_all_lead_gene,"class"]<-"MN leading gene"
sele_GO_zscore[mSC_all_lead_gene,"class"]<-"mSC leading gene"
sele_GO_zscore[intersect(MN_all_lead_gene, mSC_all_lead_gene),"class"]<-"both"
ggplot(sele_GO_zscore,aes(x=`zscore.MN`,y=`zscore.mSC`)) + geom_point(shape=19,alpha=0.8,aes(colour=as.factor(class))) +
  #scale_colour_brewer(palette = "Set1")+
  #scale_colour_gradientn(colours = terrain.colors(10)) +
  # geom_smooth(method="lm",show.legend = TRUE,level=0.95) +
  labs(title="Z score for leading genes GO_cytosolic_ribosome")+
  xlab("Z score in MN") + ylab("Z score in mSC") +
  theme(plot.title = element_text(hjust = 0.5))
