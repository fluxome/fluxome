#!/bin/bash
#Info for scrb-seq https://scg-lib-structs.readthedocs.io/en/latest/ge/mcSCRB-seq.html

#Get fastqs
fasterq-dump --split-files SRR3290185.sra
fasterq-dump --split-files SRR3290188.sra
fasterq-dump --split-files SRR3290190.sra
fasterq-dump --split-files SRR3290193.sra
fasterq-dump --split-files SRR3290196.sra
fasterq-dump --split-files SRR3290199.sra
fasterq-dump --split-files SRR3290202.sra
fasterq-dump --split-files SRR3290203.sra
fasterq-dump --split-files SRR3290204.sra

#Get whitelist for scrb-seq
wget  https://s3.amazonaws.com/pr-journal/uhyjf6e.txt
mv uhyjf6e.txt mcSCRB-seq_oligodT.txt
tail -n +2 mcSCRB-seq_oligodT.txt | cut -f 3 > mcSCRB-seq_whitelist.txt


STAR --runThreadN 8 \
     --genomeDir /home/tchari/CellDynamicsData/reference/STAR/ \
     --outFileNamePrefix /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290185_outs/ \
     --readFilesIn /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290185_2.fastq /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290185_1.fastq \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 6 --soloUMIstart 7 --soloUMIlen 10 \
     --soloCBwhitelist /home/tchari/CellDynamicsData/harissa_timecourse/mcSCRB-seq_whitelist.txt \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate \
     --soloFeatures Gene Velocyto \
     --soloBarcodeReadLength 0

STAR --runThreadN 8 \
     --genomeDir /home/tchari/CellDynamicsData/reference/STAR/ \
     --outFileNamePrefix /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290188_outs/ \
     --readFilesIn /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290188_2.fastq /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290188_1.fastq \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 6 --soloUMIstart 7 --soloUMIlen 10 \
     --soloCBwhitelist /home/tchari/CellDynamicsData/harissa_timecourse/mcSCRB-seq_whitelist.txt \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate \
     --soloFeatures Gene Velocyto \
     --soloBarcodeReadLength 0

STAR --runThreadN 8 \
     --genomeDir /home/tchari/CellDynamicsData/reference/STAR/ \
     --outFileNamePrefix /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290190_outs/ \
     --readFilesIn /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290190_2.fastq /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290190_1.fastq \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 6 --soloUMIstart 7 --soloUMIlen 10 \
     --soloCBwhitelist /home/tchari/CellDynamicsData/harissa_timecourse/mcSCRB-seq_whitelist.txt \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate \
     --soloFeatures Gene Velocyto \
     --soloBarcodeReadLength 0

STAR --runThreadN 8 \
     --genomeDir /home/tchari/CellDynamicsData/reference/STAR/ \
     --outFileNamePrefix /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290193_outs/ \
     --readFilesIn /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290193_2.fastq /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290193_1.fastq \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 6 --soloUMIstart 7 --soloUMIlen 10 \
     --soloCBwhitelist /home/tchari/CellDynamicsData/harissa_timecourse/mcSCRB-seq_whitelist.txt \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate \
     --soloFeatures Gene Velocyto \
     --soloBarcodeReadLength 0


STAR --runThreadN 8 \
     --genomeDir /home/tchari/CellDynamicsData/reference/STAR/ \
     --outFileNamePrefix /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290196_outs/ \
     --readFilesIn /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290196_2.fastq /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290196_1.fastq \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 6 --soloUMIstart 7 --soloUMIlen 10 \
     --soloCBwhitelist /home/tchari/CellDynamicsData/harissa_timecourse/mcSCRB-seq_whitelist.txt \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate \
     --soloFeatures Gene Velocyto \
     --soloBarcodeReadLength 0

STAR --runThreadN 8 \
     --genomeDir /home/tchari/CellDynamicsData/reference/STAR/ \
     --outFileNamePrefix /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290199_outs/ \
     --readFilesIn /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290199_2.fastq /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290199_1.fastq \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 6 --soloUMIstart 7 --soloUMIlen 10 \
     --soloCBwhitelist /home/tchari/CellDynamicsData/harissa_timecourse/mcSCRB-seq_whitelist.txt \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate \
     --soloFeatures Gene Velocyto \
     --soloBarcodeReadLength 0

STAR --runThreadN 8 \
     --genomeDir /home/tchari/CellDynamicsData/reference/STAR/ \
     --outFileNamePrefix /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290202_outs/ \
     --readFilesIn /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290202_2.fastq /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290202_1.fastq \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 6 --soloUMIstart 7 --soloUMIlen 10 \
     --soloCBwhitelist /home/tchari/CellDynamicsData/harissa_timecourse/mcSCRB-seq_whitelist.txt \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate \
     --soloFeatures Gene Velocyto \
     --soloBarcodeReadLength 0

STAR --runThreadN 8 \
     --genomeDir /home/tchari/CellDynamicsData/reference/STAR/ \
     --outFileNamePrefix /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290203_outs/ \
     --readFilesIn /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290203_2.fastq /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290203_1.fastq \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 6 --soloUMIstart 7 --soloUMIlen 10 \
     --soloCBwhitelist /home/tchari/CellDynamicsData/harissa_timecourse/mcSCRB-seq_whitelist.txt \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate \
     --soloFeatures Gene Velocyto \
     --soloBarcodeReadLength 0

STAR --runThreadN 8 \
     --genomeDir /home/tchari/CellDynamicsData/reference/STAR/ \
     --outFileNamePrefix /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290204_outs/ \
     --readFilesIn /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290204_2.fastq /home/tchari/CellDynamicsData/harissa_timecourse/SRR3290204_1.fastq \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 6 --soloUMIstart 7 --soloUMIlen 10 \
     --soloCBwhitelist /home/tchari/CellDynamicsData/harissa_timecourse/mcSCRB-seq_whitelist.txt \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate \
     --soloFeatures Gene Velocyto \
     --soloBarcodeReadLength 0



