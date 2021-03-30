"""
Visualization of the distribution of mean psi in different cells.
  Input: all_psi_new, psi of each AS in each cell
"""

psi_new<-all_psi_new[which(all_psi_new$cell_type=="MN"),]
stack_data<-data.frame(matrix(NA,length(psi_new[,1])*5,5))
colnames(stack_data)<-c("cell","psi_range","ratio","mean_psi_genes","glbl_mean_psi")
num=1
for (i in rownames(MN_se_psi_mean)){
  stack_data$cell[num:(num+4)] <- i
  stack_data$mean_psi_genes[num:(num+4)] <- psi_new$mean_psi_genes[which(rownames(MN_se_psi_mean)==i)]
  stack_data$glbl_mean_psi[num:(num+4)] <- psi_new$glbl_mean_psi[which(rownames(MN_se_psi_mean)==i)]
  num=num+5
}

stack_data$psi_range<-(as.numeric(rownames(stack_data))-1) %%5 +1

num=1
for (i in 1:length(MN_se_psi_mean[,1])){
  tt<-sum(!is.na(MN_se_psi_mean[i,-length(MN_se_psi_mean[1,])]))
  stack_data$ratio[num]<-length(which(MN_se_psi_mean[i,-length(MN_se_psi_mean[1,])]==0))/tt
  stack_data$ratio[num+1]<-length(which(MN_se_psi_mean[i,-length(MN_se_psi_mean[1,])]>0 & MN_se_psi_mean[i,-length(MN_se_psi_mean[1,])]<=0.33))/tt
  stack_data$ratio[num+2]<-length(which(MN_se_psi_mean[i,-length(MN_se_psi_mean[1,])]>0.33 & MN_se_psi_mean[i,-length(MN_se_psi_mean[1,])]<=0.66))/tt
  stack_data$ratio[num+3]<-length(which(MN_se_psi_mean[i,-length(MN_se_psi_mean[1,])]>0.66 & MN_se_psi_mean[i,-length(MN_se_psi_mean[1,])]<1))/tt
  stack_data$ratio[num+4]<-length(which(MN_se_psi_mean[i,-length(MN_se_psi_mean[1,])]==1))/tt
  num=num+5
  }


ggplot(data =stack_data) + 
  geom_bar( aes(x = cell, y = ratio,fill = factor(psi_range)),
                    stat = "identity",
                    position = "fill") +
  theme(axis.text.x = element_blank()) + #°ÑxÖá±êÇ©È¥µô
  geom_line(aes(x = cell, y =mean_psi_genes,group=1))+
  geom_line(aes(x = cell, y =glbl_mean_psi,group=1))+
  geom_point(aes(x = cell, y =mean_psi_genes))+
  geom_point(aes(x = cell, y =glbl_mean_psi))

ggplot() + 
  geom_line(data =small_data,aes(x = cell, y =mean_gene_psi, group=1))

layout(matrix(1:1, 1, 4,byrow=T))

