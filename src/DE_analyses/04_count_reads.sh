#!/usr/bin/env bash
# Here we will count the reads using the featureCounts tool included in the Subread package.

baseDir="$HOME/archive/projects/PRJ0000002_CDMON"
outputDir="$baseDir/output"
#bamDir="$outputDir/bam"
bamDir="$HOME/bam"

countDir="$outputDir/counts"
mkdir $countDir

annotation="$HOME/archive/common_data/genome/homo_sapiens/ensembl_v95/Homo_sapiens.GRCh38.95.gtf"

bam_files=$(ls $bamDir/*.bam) 

#echo $bam_files

featureCounts -T 8 -a $annotation -t exon -g gene_id -o $countDir/counts.txt $bam_files
