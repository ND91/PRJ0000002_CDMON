#!/usr/bin/env bash

baseDir="$HOME/archive/projects/PRJ0000002_CDMON"
samples="$baseDir/data/samples/samplesheet_gse107011_PRJ0000002_CDMON_V1.csv"

mapfile -t urls < <(grep -Po 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR[0-9]{3}/[0-9]{3}/SRR[0-9]{7}/SRR[0-9]{7}' "$samples" | sort | uniq)
for url in ${urls[@]}; do
	url_r1=$url"_1.fastq.gz"
	echo $url_r1
	wget -P /home/ayliyim/archive/projects/PRJ0000002_CDMON/data/GSE107019 $url_r1 

	url_r2=$url"_2.fastq.gz"
	echo $url_r2
	wget -P /home/ayliyim/archive/projects/PRJ0000002_CDMON/data/GSE107019 $url_r2

	#mkdir "$gse_dir"
done