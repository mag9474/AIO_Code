#!/bin/bash
#SBATCH --job-name=run_rna
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=
#SBATCH --partition=cpu_short
#SBATCH --ntasks=1
#SBATCH --mem=30gb
#SBATCH --output=logs/align_test_%j.log
#SBATCH --time 12:00:00
#SBATCH --array=1-6

module load trimgalore/0.5.0
module load fastqc/0.11.7

# module load samtools/1.9

module load igenome

####EDIT_THIS#####
files_dir=/gpfs/scratch/mag9474/aio_folder/single_end_fastq_only/

#project dir
project_dir=$(pwd)

#setting genomeDir
genomedir=/gpfs/data/igorlab/ref/hg38/STAR.gencode.v34

#############################################################
sample=$(awk "NR==${SLURM_ARRAY_TASK_ID} {print \$1}" file_names.txt)

#############################################################

#Creating subdir and file dir prep for step1
extension=.fastq.gz
file=$files_dir${sample}$extension
output_dir_1=${project_dir}/step_1_out/
mkdir -p $output_dir_1
echo hello $file

#running trim galore
trim_galore $file -o $output_dir_1 --length 15 --fastqc --small_rna --fastqc_args "-o $output_dir_1"

#Creating subdir and file dir prep for step2
output_dir_2=${project_dir}/step_2_out/
star_extension=_trimmed.fq.gz
star_file=$output_dir_1${sample}$star_extension
mkdir -p $output_dir_2${sample}

#moving to step2 folder
cd $output_dir_2${sample}

#modules
module unload trimgalore/0.5.0
module load star/2.7.3a

# running star
echo hello $star_file

STAR --genomeDir $genomedir --readFilesIn $star_file  --readFilesCommand zcat               \
    --outFilterMultimapNmax 20 --alignIntronMax 1   --outFilterMismatchNoverLmax 0.03       \
    --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 16

#copy and rename output
cp Aligned.out.sam $output_dir_2${sample}.out.sam

module unload star/2.7.3a
module load samtools/1.9

cd $project_dir
output_dir_3=${project_dir}/step_3_out/
mkdir -p $output_dir_3/raw_bam
mkdir -p $output_dir_3/q30_bam
mkdir -p $output_dir_3/q30_sorted_bam

#converting to bam
samtools view -S -b $output_dir_2${sample}.out.sam > ${output_dir_3}raw_bam/${sample}.bam

#remove reads <30
samtools view -b -q 30 ${output_dir_3}raw_bam/${sample}.bam -o ${output_dir_3}q30_bam/${sample}_q30.bam

#sort bam file
samtools sort ${output_dir_3}q30_bam/${sample}_q30.bam  -o ${output_dir_3}q30_sorted_bam/${sample}_q30_sorted.bam

#index file
samtools index -b ${output_dir_3}q30_sorted_bam/${sample}_q30_sorted.bam

# module unload python
# module load deeptools
#
# output_dir_4=${project_dir}/step_4_out/
# mkdir -p $output_dir_4
#
# bamCoverage -b ${output_dir_3}q30_sorted_bam/${sample}_q30_sorted.bam --binSize 1 --normalizeUsing RPKM -o  ${output_dir_4}${sample}_RPKM.bw
#
# #Use htseq to get gene counts for each sample
# module unload python
# module load python/cpu/2.7.15-ES
#
# output_dir_5=${project_dir}/step_5_out/
# mkdir -p $output_dir_5
#
# htseq-count -f bam -r pos -s reverse ${output_dir_3}q30_sorted_bam/${sample}_q30_sorted.bam /gpfs/share/apps/iGenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf > $output_dir_5${sample}_counts.txt

#
