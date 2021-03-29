#!/bin/bash -l
#SBATCH -J job_atacpipeline
#SBATCH -N 1
#SBATCH --partition=fn_medium


module load bowtie2/2.3.4.1
module load samtools/1.9
module load bedtools
module load igvtools/2.4.1
module load igenome

cd /gpfs/data/courses/bmscga4498/2021/TE_analysis

#############################################################
sample=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" file.txt)

bowtie2 --no-discordant -p 12 --no-mixed -N 0  --un unaligned_${sample}.sam -x /gpfs/home/sb5169/sb5169/hg38/hg38 -1 /gpfs/data/courses/bmscga4498/2021/RabbaniRun/FASTQ-CLEAN/${sample}_R1.fastq.gz -2 /gpfs/data/courses/bmscga4498/2021/RabbaniRun/FASTQ-CLEAN/${sample}_R2.fastq.gz -S ${sample}.sam > bowtie_${sample}.outout

#convert sam to bam
samtools view -b  -h ${sample}.sam -o ${sample}.bam

#sort bam
samtools sort ${sample}.bam -o ${sample}_sorted.bam

###Index BAM file to obtain .bai files
samtools index -b ${sample}_sorted.bam

#counts
bedtools multicov -bams ${sample}_sorted.bam -bed hg38_removeclasses.bed > ${sample}_bowtie_elements_counts.txt
awk '{print $1":"$2":"$3":"$4"\t"$5}' ${sample}_bowtie_elements_counts.txt   > ${sample}_bowtie_elements_counts_2.txt  


#samtools view -h -o ${sample}_sorted.sam ${sample}_sorted.bam
#igvtools count -w 50 -e 250 ${sample}_sorted.bam ${sample}_sorted.tdf /gpfs/home/sb5169/sb5169/hg38/hg38.fa


