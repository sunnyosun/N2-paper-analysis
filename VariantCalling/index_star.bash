#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#module load bwa
#bwa index /ifs/data/proteomics/projects/Sunny/genome/hg38.fa

#module load bowtie2
#genome=$1
#bowtie2-build /ifs/data/proteomics/projects/Sunny/rnaseq/rnaseq_genome/${genome}.fa ${genome}

# 1-.fa
# 2-.gtf
# 3- index directory

module load star/2.4.5a

STAR --runMode genomeGenerate --runThreadN 8 --genomeDir /ifs/data/proteomics/projects/Sunny/jingchuan/yeast_genome/BY4741_star --genomeFastaFiles BY4741_genome.fa --sjdbGTFfile BY4741_genome.gff --outFileNamePrefix BY4741_genome --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS --genomeSAindexNbases 10

