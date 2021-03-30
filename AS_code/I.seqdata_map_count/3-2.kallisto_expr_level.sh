#!/bin/sh

export PATH=/home/yihanlin_pkuhpc/lustre3/yihanlin_cls/cxy/cxyapp/kallisto:$PATH
source /appsnew/source/bioapps/sratoolkit.2.10.sh

seq_list=`ls ./raw2 | grep "SRR"`

input_dir="./raw2/"
output_dir="./fastq/"
#mkdir ${output_dir}
#for i in ${seq_list};
#do
#       fasterq-dump --split-3 ${input_dir}${i} 
#done

input_dir="./fastq/"
input_fasta="./genome_transcript/gencode.v36.transcripts.fa"
input_gtf="./genome_transcript/gencode.v36.annotation.gtf"
ERCC_fasta="./genome_transcript/gencode.v36.transcripts.ERCC.fa"
ERCC_gtf="./genome_transcript/gencode.v36.annotation.ERCC.gtf"
cp ${input_fasta} ${ERCC_fasta}
cp ${input_gtf} ${ERCC_gtf}
cat "./ERCC92/ERCC92.fa" >> ${ERCC_fasta}
cat "./ERCC92/ERCC92.gtf" >> ${ERCC_gtf}
output_dir="./kal_quant/"
index_dir="kal_index"

kallisto index -i ${index_dir} ${ERCC_fasta}
for i in ${seq_list};
do
	kallisto quant -i ${index_dir} -o ${output_dir}${i} -g ${ERCC_gtf} ${input_dir}${i}_1.fastq  ${input_dir}${i}_2.fastq
done

python kal_sc_aggregate.py
