#!/usr/bin/env bash

dir_name="${1:-PRJNA1066021}"

mkdir -p ${dir_name}

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/028/SRR27599428/SRR27599428_1.fastq.gz -O ${dir_name}/PAO1_control_1__1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/021/SRR27599421/SRR27599421_1.fastq.gz -O ${dir_name}/PAO1_MEM_3__1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/024/SRR27599424/SRR27599424_1.fastq.gz -O ${dir_name}/PAO1_control_3__1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/022/SRR27599422/SRR27599422_1.fastq.gz -O ${dir_name}/PAO1_MEM_2__1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/027/SRR27599427/SRR27599427_1.fastq.gz -O ${dir_name}/PAO1_control_2__1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/028/SRR27599428/SRR27599428_2.fastq.gz -O ${dir_name}/PAO1_control_1__2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/022/SRR27599422/SRR27599422_2.fastq.gz -O ${dir_name}/PAO1_MEM_2__2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/021/SRR27599421/SRR27599421_2.fastq.gz -O ${dir_name}/PAO1_MEM_3__2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/024/SRR27599424/SRR27599424_2.fastq.gz -O ${dir_name}/PAO1_control_3__2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/027/SRR27599427/SRR27599427_2.fastq.gz -O ${dir_name}/PAO1_control_2__2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/023/SRR27599423/SRR27599423_1.fastq.gz -O ${dir_name}/PAO1_MEM_1__1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR275/023/SRR27599423/SRR27599423_2.fastq.gz -O ${dir_name}/PAO1_MEM_1__2.fastq.gz
