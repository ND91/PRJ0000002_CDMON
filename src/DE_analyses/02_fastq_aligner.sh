#!/usr/bin/env bash
# Aligns the fastq.gz files against the GRCh38 Ensembl genome annotated using v95

baseDir="$HOME/archive/projects/PRJ0000002_CDMON"
fastqDir="$HOME/archive/raw_data/fastq/rnaseq/NCBI_GSE107019"
outputDir="$baseDir/output"

genomeDir="$HOME/ensembl_GRCh37_r87_sjdbO49"

samDir="$HOME/sam"
mkdir $samDir

samplesheet="$baseDir/data/samples/samplesheet_gse107011_PRJ0000002_CDMON_V1.csv"
srrs=($(gawk 'match($0, /(SRR[0-9]+),.+$/, arr) {print arr[1]}' $samplesheet))

STAR --genomeDir $genomeDir --genomeLoad LoadAndExit |& tee -a $samDir/STAR.log.txt

for srr in "${srrs[@]}"; do
	echo "Starting on $srr"
	    
	sampleDir="$samDir/$srr"
	mkdir $sampleDir
	fastq_1="$fastqDir/${srr}_1.fastq.gz"
	fastq_2="$fastqDir/${srr}_2.fastq.gz"
	
	if [ -f $fastq_1 ]; then
		echo "$srr R1 found"
    else
		echo "$srr R1 not found"
		continue
	fi
	
	if [ -f $fastq_2 ]; then
		echo "$srr R2 found"
    else
		echo "$srr R2 not found"
		continue
	fi
	
	STAR --runThreadN 8 --genomeDir $genomeDir --genomeLoad LoadAndKeep --readFilesIn $fastq_1 $fastq_2 --readFilesCommand zcat --outFileNamePrefix $sampleDir/$srr_ |& tee -a $samDir/STAR.log.txt
done

STAR --genomeDir $genomeDir --genomeLoad Remove |& tee -a $samDir/STAR.log.txt
