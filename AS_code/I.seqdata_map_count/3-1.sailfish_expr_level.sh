#!/bin/sh

#sailfish
export PATH=/home/yihanlin_pkuhpc/lustre3/yihanlin_cls/cxy/cxyapp/SailfishBeta-0.10.0_CentOS5/bin:$PATH

prefix="./iPSC"
seq_list=`cat ${prefix}list`

outdir="${prefix}_filted_fa/"
index="./trans_index"
mkdir ${outdir}

#build index ran once
#sailfish index -t ./genome/gencode.v19.pc_transcripts.fa -o ${index}
for i in ${seq_list};
do
        #if  [ -e ${output_dir}SRR4047${i}_R1.fa ]; then
        #        echo "cut SRR4047${i}"
        #else
        #     
	samtools fasta -F 0x900 -1 ${outdir}SRR4047${i}R1.fa -2 ${outdir}SRR4047${i}R2.fa ${prefix}_filt_mapped/SRR4047${i}Aligned.out.bam
	sailfish quant -p 20 -i ${index} -l IU -1 ${outdir}SRR4047${i}R1.fa -2 ${outdir}SRR4047${i}R2.fa -o ${outdir}SRR4047${i}
# -p number of threads 
	#fi
done

rm ${outdir}*.fa

