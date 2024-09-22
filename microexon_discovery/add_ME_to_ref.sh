#!/bin/bash

#run script within the module_3 conda env
#firstly convert extended.reference.gtf file to bed12 format

filtered_ME=$1 #"Report/out.high_quality.txt"
ext_gtf=$2 #'/g/korbel/zafferani/RNA_seq_snakepipes/integrated_patients/gffcompare/extended_ref_annotation.gtf'
Genome_fasta=$3 #'/g/korbel/zafferani/genome_indexed/human_gn/Homo_sapiens.GRCh38.dna.primary_assembly.fa'


bed12='Report/extended_ref_annotation.bed12'

bedparse gtf2bed $ext_gtf > $bed12 && wait

chrM='False'

ME_GTF="Report/out.high_quality.gtf"

#convert microexons text file to GTF format and append them to extended reference GTF
python get_isoforms2.py $Genome_fasta $bed12 $ext_gtf $filtered_ME $chrM  > $ME_GTF && wait

if [ -e "$ME_GTF" ]; then
    if [ -s "$ME_GTF" ]; then
        echo "File "$ME_GTF" exists and is not empty."
    else
        echo "File "$ME_GTF" exists but is empty."
        rm "$ME_GTF"
        exit -1
    fi
else
    echo "File "$ME_GTF" does not exist."
    exit -1
fi

final_output='Report/out.high_quality_calibrated.gtf'

#filtering the transcriptome extended with microexons
python final_filtering_GTF.py "$ME_GTF" "$final_output" && wait

if [ -e "$final_output" ]; then
    if [ -s "$final_output" ]; then
        echo "File "$final_output" exists and is not empty."
    else
        echo "File "$final_output" exists but is empty."
        rm  "$final_output"
        exit -1
    fi
else
    echo "File "$final_output" does not exist."
    exit -1
fi