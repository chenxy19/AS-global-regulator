#!/bin/sh

source /appsnew/source/bioapps/sratoolkit.2.10.sh
#fastQC
export PATH=/home/yihanlin_pkuhpc/lustre2/cxy/apps/fastqc/FastQC:/lustre2/yihanlin_pkuhpc/cxy/apps/regtools/build:$PATH
#cutadapt
export PATH=/home/yihanlin_pkuhpc/yihanlin_cls/.local/bin:$PATH
#STAR
export PATH=/home/yihanlin_pkuhpc/lustre3/yihanlin_cls/cxy/cxyapp/STAR-2.7.5a/source:$PATH
#sailfish
#export PATH=/home/yihanlin_pkuhpc/lustre3/yihanlin_cls/cxy/cxyapp/SailfishBeta-0.10.0_CentOS5/bin:$PATH
#Seqmonk
export PATH=/home/yihanlin_pkuhpc/lustre3/yihanlin_cls/cxy/cxyapp/SeqMonk:$PATH
#MISO
#export PATH=/home/yihanlin_pkuhpc/lustre3/yihanlin_cls/cxy/AS_data/miso-env/bin:$PATH


seq_list=`ls ./raw | grep "ERR"`

input_dir="./raw/"
output_dir="./fastq/"
mkdir ${output_dir}
for i in ${seq_list};
do
       fasterq-dump --threads ${nthread} --split-3 ${input_dir}${i} -O ${output_dir}
done

################## FastQC #####################

FastQC_outdir="./FastQC_pre_trim/"
mkdir "${FastQC_outdir}"
FASTQ_dir="./fastq/"
ls -1 "${FASTQ_dir}" > FASTQ_files
#FASTQ_files=`cat FASTQ_files`
FASTQ_files=`head FASTQ_files`
for i in ${FASTQ_files};
do
       echo "############################################"
       date
       echo "FastQC processing file: ${FASTQ_dir}${i}"
       echo "############################################"
       fastqc --outdir "${FastQC_outdir}" --threads 12 "${FASTQ_dir}${i}"
       echo
done

#FastQC_outdir="./FastQC_post_trim/"
#mkdir "${FastQC_outdir}"
#FASTQ_dir="./cut/"
#ls -1 "${FASTQ_dir}" > FASTQ_files
#FASTQ_files=`cat FASTQ_files`
#for i in ${FASTQ_files};
#do
#       echo "############################################"
#       date
#       echo "FastQC processing file: ${FASTQ_dir}${i}"
#       echo "############################################"
#       fastqc --outdir "${FastQC_outdir}" --threads 12 "${FASTQ_dir}${i}"
#       echo
#done
