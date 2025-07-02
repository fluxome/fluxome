#!/bin/bash
#Build star index file for mouse mm10

mkdir -p /home/tchari/CellDynamicsData/reference/STAR
STAR --runThreadN 20 --runMode genomeGenerate --genomeSAindexNbases 13 --genomeDir /home/tchari/CellDynamicsData/reference/STAR --genomeFastaFiles /home/tchari/CellDynamicsData/reference/mm10.fa --sjdbGTFfile /home/tchari/CellDynamicsData/reference/mm10.refGene.gtf
