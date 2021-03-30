import pandas as pd
import os
import re

path = "kal_quant/" #文件夹目录
files= os.listdir(path)

i=0
for folder in files:
    try:
        aff=pd.read_csv(path+folder+"/abundance.tsv",sep="\t",header=0,index_col=0)
        aff[f"tpm_{folder}"]=aff["tpm"]
        aff[f"est_counts_{folder}"]=aff["est_counts"]
        if i==0:
            TPM=aff[f"tpm_{folder}"].to_frame()
            nreads= aff[f"est_counts_{folder}"].to_frame()
        else:
            TPM=pd.merge(TPM, aff[f"tpm_{folder}"].to_frame(),how='inner',left_index=True,right_index=True)
            nreads=pd.merge(nreads, aff[f"est_counts_{folder}"].to_frame(),how='inner',left_index=True,right_index=True)
        i=i+1
        
    except:
        pass

ENSG=[]
for name in aff.index:
    try:
        ensg=re.search("^(ENST.*?)\|(ENSG.*?)\|",name).group(2)
        ENSG.append(ensg)
    except:
        ENSG.append(name)
        print(name)
TPM.index=ENSG
nreads.index=ENSG

TPM.to_csv(path+"hHSC_kal_TPM.csv")
nreads.to_csv(path+"hHSC_kal_nreads.csv")
