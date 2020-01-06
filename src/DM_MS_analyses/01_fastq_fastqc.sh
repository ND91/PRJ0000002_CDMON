#!/usr/bin/env bash
# Run FastQC on the .fastq files to ensure the files are of good quality

baseDir="$HOME/archive/projects/PRJ0000002_CDMON"
fastqDir="$HOME/archive/raw_data/fastq/methylseq/NCBI_GSE73788"

fastqcDir="$baseDir/output/fastqc"
mkdir $fastqcDir

samplesheet="$baseDir/data/samples/samplesheet_methylseq_PRJ0000002_CDMON_V1.csv"

srrs=($(gawk 'match($0, /(SRR[0-9]+),.+$/, arr) {print arr[1]}' $samplesheet))

for srr in "${srrs[@]}"; do
	#echo $srr
	
 	sampleDir="$fastqcDir/$srr"
 	mkdir $sampleDir

 	fastq="$fastqDir/$srr.fastq.gz"
	#echo $fastq
	if [ -f $fastq ]; then
		echo "$srr found"
	else
		"$srr not found"
		continue
	fi

	fastqc --outdir="$sampleDir" --threads 8 $fastq |& tee -a "$fastqcDir/fastqc.log.txt"

done
