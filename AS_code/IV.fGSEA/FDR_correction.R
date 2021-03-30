
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "FDR_P_out.csv"
}

p_value=read.csv(args[1],stringsAsFactors = F)
#p_value$BH_p<-p.adjust(p_value$p, "BH")
p_value$BH_p<-p.adjust(p_value[,2], "BH")
#print(p_value[,1])
write.table (p_value[which(p_value$BH_p<0.05),], file = args[2], row.names =TRUE, col.names =TRUE, quote =TRUE)

