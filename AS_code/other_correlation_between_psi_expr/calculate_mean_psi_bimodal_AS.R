"""
calculate 2 types of mean psi only in bimodal genes (the version where each run is a sample)
  Input: se_bi_psi.csv (colnames are bimodal events), events.csv, reads.csv
  Output: sum_junc into junc_sum_bimodal.csv (3 rows * n_sample), 
          sam_psi into junc_psi_bimodal.csv (AS events * n_sample),
          iPSC_bi_2psi (n_sample * 2rows: bi_glbl_mean_psi and bi_gene_mean_psi) 
          
"""
########## 1. calculate global psi for bimodal genes ###############
# only for bimodal genes
#"mouseSC/mouseSC_",
for (prefix in c("/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/iPSC_")){
  
  bimodal_file = paste(prefix,"se_bi_psi.csv",sep="")
  bimodal= read.csv(bimodal_file,header = TRUE,stringsAsFactors = FALSE) 
  
  event_file=paste(prefix,"mean_psi/events.csv",sep="")
  events= read.csv(event_file,header = TRUE,stringsAsFactors = FALSE) 
  junctions<-events[,c("event_id","junction13","junction12","junction23")]    #其实unique只有8690,which is reasonable
  junctions<-unique(junctions)
  junctions$event_id<- gsub("[^0-9A-z]",".",junctions$event_id)
  junctions <- junctions[which(junctions$event_id %in% colnames(bimodal)),]
  
  junc_reads_file=paste(prefix,"mean_psi/reads.csv",sep="")
  junc_reads= read.csv(junc_reads_file,header = TRUE,stringsAsFactors = FALSE) 
  head(junc_reads)
  
  SRRlist<-unique(junc_reads$sample_id)
  sum_junc<-data.frame(matrix(NA,3,length(SRRlist)))
  colnames(sum_junc)<-SRRlist
  rownames(sum_junc)<-c("junc13","junc12","junc23")
  ASlist<-junctions$event_id
  
  sam_psi<-data.frame(matrix(NA,nrow(junctions),length(SRRlist)))
  colnames(sam_psi)<-SRRlist
  rownames(sam_psi)<-junctions$event_id
  
  #对于SRR1
  for (j in 1:length(SRRlist)){
    junctions$read13 <- 0
    junctions$read12 <- 0
    junctions$read23 <- 0
    cell_junc<-junc_reads[which(junc_reads$sample_id==SRRlist[j]),]
    for (i in 1:length(junctions[,1])){   #不管正负链，永远是小的在前，大的在后的
      #junc_reads$sample_id[which(junc_reads$junction_id==junctions$junction13[i])]
      #不需要每个junction必须有reads!这里太苛刻了！
      if (length(which(cell_junc$junction_id==junctions$junction13[i]))!=0){
        junctions$read13[i]<-cell_junc$reads[which(cell_junc$junction_id==junctions$junction13[i])]
      }
      if (length(which(cell_junc$junction_id==junctions$junction12[i]))!=0){
        junctions$read12[i]<-cell_junc$reads[which(cell_junc$junction_id==junctions$junction12[i])]
      }
      if (length(which(cell_junc$junction_id==junctions$junction23[i]))!=0){
        junctions$read23[i]<-cell_junc$reads[which(cell_junc$junction_id==junctions$junction23[i])]
      }
    }
    adeq_junc<-which(junctions$read13 + junctions$read12 + junctions$read23 >= 10)  #从134个升到169，提升并不显著
    sum_junc[1,j]<-sum(junctions[adeq_junc,"read13"])
    sum_junc[2,j]<-sum(junctions[adeq_junc,"read12"])
    sum_junc[3,j]<-sum(junctions[adeq_junc,"read23"])
    
    for (k in adeq_junc){
      sam_psi[k,SRRlist[j]]<-(junctions$read12[k] + junctions$read23[k])/(junctions$read12[k] + junctions$read23[k] + 2*junctions$read13[k])
    }
    print(SRRlist[j])
  }
  junc_mat_file=paste(prefix,"junc_sum_bimodal.csv",sep="")
  write.table (sum_junc, file = junc_mat_file, row.names =TRUE, col.names =TRUE, quote =TRUE)
  junc_mat_file=paste(prefix,"junc_psi_bimodal.csv",sep="")
  write.table (sam_psi, file = junc_mat_file, row.names =TRUE, col.names =TRUE, quote =TRUE)
  try_globl_mean_psi<-as.data.frame(t(sum_junc))
  try_globl_mean_psi$mean_psi<-(try_globl_mean_psi$junc12 + try_globl_mean_psi$junc23)/(try_globl_mean_psi$junc12 + try_globl_mean_psi$junc23 + 2*try_globl_mean_psi$junc13)

  try_sam_psi<-as.data.frame(t(sam_psi))
  try_sam_psi$gene_mean_psi <- apply(try_sam_psi, 1, mean_na_rm)
  
  iPSC_bi_2psi <- merge_by_rownames(try_globl_mean_psi, try_sam_psi)
  iPSC_bi_2psi<-iPSC_bi_2psi[,c("mean_psi","gene_mean_psi")]
  colnames(iPSC_bi_2psi)<-c("bi_glbl_mean_psi","bi_gene_mean_psi")
  junc_mat_file=paste(prefix,"psi_bimodal.csv",sep="")
  write.table (iPSC_bi_2psi, file = junc_mat_file, row.names =TRUE, col.names =TRUE, quote =TRUE)
  
  
}


iPSC2p_2meth<-merge_by_rownames(iPSC_bi_2psi,as.data.frame(iPSC_glbl_mean_psi))
iPSC2p_2meth<-merge_by_rownames(iPSC2p_2meth, as.data.frame(iPSC_gene_mean_psi))

write.table (iPSC2p_2meth, file = "/Users/chenxinyi/Desktop/output/data1_MN+NPC+iPSC/iPSC_2p_2meth.csv", row.names =TRUE, col.names =TRUE, quote =TRUE)

#### 经过验证，用旧的代码，或者上述函数，求全部基因的psi能得到一样的结果～

mean_na_rm<-function(a){
  return(mean(a, na.rm = TRUE))
}


# 确认一下这两段代码除此之外没有其他区别！
plot(iPSC2p_2meth$bi_glbl_mean_psi,iPSC2p_2meth$bi_gene_mean_psi,main="iPSC bimodal AS",xlab="mean psi 1",ylab="mean psi 2")
LM<-lm(iPSC2p_2meth$iPSC_glbl_mean_psi~iPSC2p_2meth$iPSC_gene_mean_psi)
plot(iPSC2p_2meth$iPSC_glbl_mean_psi,iPSC2p_2meth$iPSC_gene_mean_psi, main="iPSC all AS",xlab="mean psi 1",ylab="mean psi 2")

