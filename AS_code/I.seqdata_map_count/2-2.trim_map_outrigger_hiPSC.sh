#!/bin/sh

#cutadapt
export PATH=/home/yihanlin_pkuhpc/yihanlin_cls/.local/bin:$PATH
#STAR
export PATH=/home/yihanlin_pkuhpc/lustre3/yihanlin_cls/cxy/cxyapp/STAR-2.7.5a/source:$PATH
#sailfish
#export PATH=/home/yihanlin_pkuhpc/lustre3/yihanlin_cls/cxy/cxyapp/SailfishBeta-0.10.0_CentOS5/bin:$PATH
#MISO
export PATH=/home/yihanlin_pkuhpc/lustre3/yihanlin_cls/cxy/AS_data/miso-env/bin:$PATH


#self-defining part:
#seq_list=$(seq 279 290)
prefix="./iPSC"
seq_list=`cat ${prefix}list`


############### Cutadapt trim the adapters #####################

input_dir="${prefix}fastq/"
output_dir="${prefix}_cut/"
mkdir ${output_dir}

for i in ${seq_list};
do
	if  [ -e ${output_dir}SRR4047${i}_R1.fa ]; then
		echo "cut SRR4047${i}"
	else
	#	 cutadapt support pair-end trimming
		cutadapt -j 0 -b TCGTATGCCGTCTTCTGCTTG -b ATCTCGTATGCCGTCTTCTGCTTG -b CGACAGGTTCAGAGTTCTACAGTCCGACGATC -b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -b "A{50}" -b "T{50}" -B TCGTATGCCGTCTTCTGCTTG -B ATCTCGTATGCCGTCTTCTGCTTG -B CGACAGGTTCAGAGTTCTACAGTCCGACGATC -B GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -B "A{50}" -B "T{50}"  -o ${output_dir}SRR4047${i}_R1.fa -p ${output_dir}SRR4047${i}_R2.fa ${input_dir}SRR4047${i}_1.fastq ${input_dir}SRR4047${i}_2.fastq  
	fi
done

################# Mask the reads to repetitive elements #########################
#build genome index for Repbase sequences
#input_fasta="./RMRBSeqs.fasta"
#output_dir="./rep_index/"
#mkdir ${output_dir}
#STAR --runThreadN 6 --limitGenomeGenerateRAM 31434910649 --runMode genomeGenerate --genomeDir ${output_dir} --genomeFastaFiles ${input_fasta}
#if limitGenomeGenerateRAM, the maximal size of genome is not set, the default threshold will be too small. The genome index only needs to be built once.


#mapping to repetitive elements
input_ref="./rep_index/"
input_dir="${prefix}_cut/"
output_dir="${prefix}_rep_filtered/"
mkdir ${output_dir}
for i in ${seq_list};
do
	if  [ -e ${output_dir}SRR4047${i}Aligned.out.sam ]; then
                echo "repetitive elements filtered SRR4047${i}"
        else        
		STAR --runThreadN 12 --genomeDir ${input_ref} --outFilterMultimapNmax 100000 --readFilesIn ${input_dir}SRR4047${i}_R1.fa ${input_dir}SRR4047${i}_R2.fa --outFileNamePrefix ${output_dir}SRR4047${i} --outReadsUnmapped Fastx #${output_dir}SRR4047${i} 
	fi
done

#Only when a read is mapped to more than 100000 loci can it be defined as unmapped. In the way, all reads mapping to repetitive elements are mapped and discarded in this step.(default is 10, while it is set to a lower level in subsequent mapping)


################# STAR mapping #####################
#build genome index for whole human genome (everything on the chromosome)
input_fasta="./genome/GRCh37.primary_assembly.genome.fa"
input_gtf="./genome/gencode.v34lift37.annotation.gtf"
output_dir="./genome_index/"
#mkdir ${output_dir}
#STAR --runThreadN 6 --runMode genomeGenerate --genomeDir ${output_dir} --genomeFastaFiles ${input_fasta} --sjdbGTFfile ${input_gtf} --sjdbOverhang 99

#input_dir="./cut_MNfastq/"
# RepeatMasker is too slow for RNA seq mapping, we will adopt its library but instead use STAR to map
#RepeatMasker -qq -pa 10 -species human ${input_dir}SRR4047279_R1.fa

#mapping to human genome
input_ref="./genome_index/"
input_dir="${prefix}_rep_filtered/"
output_dir="${prefix}_filt_mapped/"
mkdir ${output_dir}
for i in ${seq_list};
do
	if  [ -e ${output_dir}SRR4047${i}Aligned.out.sam ]; then
                echo "STAR mapped to genome SRR4047${i}"
        else
		STAR --runThreadN 12 --genomeDir ${input_ref} --outFilterMultimapNmax 5 --readFilesIn ${input_dir}SRR4047${i}Unmapped.out.mate1 ${input_dir}SRR4047${i}Unmapped.out.mate2 --outFileNamePrefix ${output_dir}SRR4047${i}
	fi
done

#multimap has been set to a very low threshold (default 10)

################# Outrigger identify the junctions and calculate psi #######################
#for outrigger
#source activate outrigger-env
out_dir="${prefix}_outrigger_5/"
mkdir ${out_dir}
outrigger index --min-reads 5 --sj-out-tab ${prefix}_filt_mapped/SRR4047*SJ.out.tab  --gtf ./genome/gencode.v34lift37.annotation.gtf --output ${out_dir}
outrigger validate --genome hg19 --fasta ./genome/GRCh37.primary_assembly.genome.fa --index ${out_dir}index -o ${out_dir}
outrigger psi --index ${out_dir}index --sj-out-tab ${prefix}_filt_mapped/SRR4047*SJ.out.tab  --ignore-multimapping -o ${out_dir}

#conda deactivate
echo "All completed for iPSC (outrigger psi)"


