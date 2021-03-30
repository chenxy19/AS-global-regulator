#!/bin/sh


source /appsnew/source/bioapps/sratoolkit.2.10.sh
#export PATH=~/.aspera/connect/bin:$PATH


#list="./SRR_List"
#ascp -v -T -l 500m -P33001 -k 1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --mode recv --host  fasp.sra.ebi.ac.uk --user era-fasp --file-list list ./fastq/
                       #-v是verbose模式，Q自适应流量控制，T取消加密，否则有的数据下载不了，i提供私钥的地址；l设置最大传输速度，不设置传输反而较慢，开高了也不好；k断点续传
#mkdir iPSCfastq/

for i in $(seq 280 347);
do
	
	#prefetch SRR4047${i} -O ./
	wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR4047${i}/SRR4047${i}.1
	#wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR81800${i}/SRR81800${i}.1
	#mv SRR81800${i}/SRR81800${i}.sra ./
        fasterq-dump --split-3 SRR4047${i}.1 -o MNfastq/SRR4047${i}
done

for i in $(seq 348 353);
do

        #prefetch SRR4047${i} -O ./
        wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR4047${i}/SRR4047${i}.1
        #wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR81800${i}/SRR81800${i}.1
        #mv SRR81800${i}/SRR81800${i}.sra ./
        fasterq-dump --split-3 SRR4047${i}.1 -o MNfastq/SRR4047${i}
done
#for i in $(seq 826 892);
#do
	#if [ -e SRR7791${i}.1 ]; then
        #        echo "acquired SRR7791${i}"
	#else                                  #防止因文件头的原因而403
        	#prefetch SRX4646${i} -O ./
#		wget -O 'test.zip' "https://sra-downloadb.st-va.ncbi.nlm.nih.gov/sos1/sra-pub-run-16/SRR7791${i}/SRR7791${i}.1"
	#mv SRR8180${i}/SRR8180${i}.sra ./
        #rm -r SRR8180${i}
       # 	fasterq-dump --split-3 SRR7791${i}.1
	#fi
#done
