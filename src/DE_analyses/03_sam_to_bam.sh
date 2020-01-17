#!/usr/bin/env bash
# Sorts and removes mitochondrial sequences as well as converts the sam to bam files

baseDir="$HOME/archive/projects/PRJ0000002_CDMON"
outputDir="$baseDir/output"
samDir="$HOME/sam"

bamDir="$HOME/bam"
mkdir $bamDir

samplesheet="$baseDir/data/samples/samplesheet_gse107011_PRJ0000002_CDMON_V1.csv"
srrs=($(gawk 'match($0, /(SRR[0-9]+),.+$/, arr) {print arr[1]}' $samplesheet))

for srr in "${srrs[@]}"; do
	echo "Starting on $srr"

	sam="$samDir/$srr/Aligned.out.sam"
	if [ -f $sam ]; then
		echo "$srr found at $sam"
	else
		echo "$srr not found at $sam"
		continue	
	fi

	#Remove unmapped reads and multiple mapings (4 + 256 = 260)
	#Remove reads with mapping score < 10
	#Remove mitochondrial sequences
	#Convert SAM to BAM
	samtools view -@ 8 -S -h -F 260 -q 10 $sam | awk '($1 ~ /^@/) || ($3 != "MT") { print $0 }' | samtools view -@ 8 -bS -o $bamDir/$srr.tmp.bam - |& tee -a $bamDir/samtools.log.txt
	#Sort reads
	samtools sort -@ 8 -m 4G -o $bamDir/$srr.bam $bamDir/$srr.tmp.bam |& tee -a $bamDir/samtools.log.txt
	#Create index
	samtools index $bamDir/$srr.bam |& tee -a $bamDir/samtools.log.txt
	#Clean up
	rm $bamDir/$srr.tmp.bam
done
