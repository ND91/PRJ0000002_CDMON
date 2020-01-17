#!/usr/bin/env bash
# Run FastQC on the .fastq files to ensure the files are of good quality

baseDir="$HOME/archive/projects/PRJ0000002_CDMON"
fastqDir="$HOME/archive/raw_data/fastq/rnaseq/NCBI_GSE107019"

fastqcDir="$baseDir/output/02_DE_analyses/fastqc"
mkdir $fastqcDir

samplesheet="$baseDir/data/samples/samplesheet_gse107011_PRJ0000002_CDMON_V1.csv"

srrs=($(gawk 'match($0, /(SRR[0-9]+),.+$/, arr) {print arr[1]}' $samplesheet))

for srr in "${srrs[@]}"; do
	#echo $srr
	
 	sampleDir="$fastqcDir/$srr"
 	mkdir $sampleDir

 	#Read 1
 	fastq_r1="$fastqDir/${srr}_1.fastq.gz"
	echo $fastq_r1
	if [ -f $fastq_r1 ]; then
		echo "$srr read 1 found"
	else
		echo "$srr read 1 not found"
		continue
	fi

	fastqc --outdir="$sampleDir" --threads 8 $fastq_r1 |& tee -a "$fastqcDir/fastqc.log.txt"

	#Read 2
	fastq_r2="$fastqDir/${srr}_2.fastq.gz"
	echo $fastq_r2
	if [ -f $fastq_r2 ]; then
		echo "$srr read 2 found"
	else
		echo "$srr read 2 not found"
		continue
	fi

	fastqc --outdir="$sampleDir" --threads 8 $fastq_r2 |& tee -a "$fastqcDir/fastqc.log.txt"

done
