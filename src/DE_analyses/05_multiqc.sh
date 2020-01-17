#!/usr/bin/env bash
# Summarize the output statistics

baseDir="$HOME/archive/projects/PRJ0000035_COLSPLSTIM"

resultsDir="$baseDir/output"

outputDir="$resultsDir/multiqc"
mkdir $outputDir

multiqc -o $outputDir $resultsDir
