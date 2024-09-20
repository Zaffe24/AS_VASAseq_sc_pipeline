#!/bin/bash


if [ $# -ne 7 ]
then
    echo "Please, give:"
    echo "(1) reference fasta file (mouse, human or full path)"
    echo "(2) fastq file to map"
    echo "(3) prefix output"
    echo "(4) Path to bwa software"
    echo "(5) Path to samtools"
    echo "(6) Stranded protocol? (y/n)"
    echo "(7) Path to riboread-selection.py script"
    exit
fi
    
ref=$1 #rRNA_human_mouse.fa
fq=$2 #/scratch/zafferani/P27_github/trimmed/EMB-JL-v102_NSXA240327_010_cbc_trimmed_homoATCG.fq.gz
out=$3 #/scratch/zafferani/P27_github/ribo_depleted/EMB-JL-v102_NSXA240327_010_cbc_trimmed_homoATCG.fq.gz
p2bwa=$4 #/g/easybuild/x86_64/Rocky/8/haswell/software/BWA/0.7.17-GCCcore-11.2.0/bin
p2samtools=$5 #/g/easybuild/x86_64/Rocky/8/haswell/software/SAMtools/1.14-GCC-11.2.0/bin
stranded=$6 #"y" 
p2s=$7 # "scripts"

folder=${fq%/*} #/scratch/zafferani/P27_github/trimmed
fq=$(basename "$fq") #EMB-JL-v102_NSXA240327_010_cbc_trimmed_homoATCG.fq.gz
out=${out%_cbc*z} #/scratch/zafferani/P27_github/ribo_depleted/EMB-JL-v102_NSXA240327_010

echo ref: $ref
echo folder: $folder 
echo fq: $fq 
echo out: $out

#checks that that fasta file has been indexed
if [ ! -f $ref.sa ]; then
    echo "Building index for reference fasta..."
    ${p2bwa}/bwa index ${ref}
fi


${p2bwa}/bwa aln ${ref} ${folder}/${fq} > ${folder}/aln_${fq%_cbc*z}.sai 
${p2bwa}/bwa samse ${ref} ${folder}/aln_${fq%_cbc*z}.sai ${folder}/${fq} | ${p2samtools}/samtools view -Sb > ${out}.aln-ribo.bam &

echo "reach this point ?"

# mapping normal reads
${p2bwa}/bwa mem -t 2 -h 15 ${ref} ${folder}/${fq} | ${p2samtools}/samtools view -Sb > ${out}.mem-ribo.bam & 

wait

${p2samtools}/samtools merge -n -r -h ${out}.aln-ribo.bam --threads 2 ${out}.all-ribo.bam ${out}.aln-ribo.bam ${out}.mem-ribo.bam 
rm ${out}.aln-ribo.bam ${out}.mem-ribo.bam ${folder}/aln_${fq%_cbc*z}.sai
${p2samtools}/samtools sort -n --threads 2 ${out}.all-ribo.bam -O BAM -o ${out}.nsorted.all-ribo.bam
rm ${out}.all-ribo.bam

python ${p2s}/riboread-selection.py ${out}.nsorted.all-ribo.bam $stranded ${out}

#exit

#if [ $stranded == "n" ]
#then
#    # select reads that don't map and create a new fastq file
#    ${p2samtools}/samtools view -f 4 ${out}.ribo.sam | awk 'BEGIN {OFS="\n"} {print "@"$1, $10, "+", $11}' > ${out}.nonRibo.fastq & 
#    # select reads that map (even low qualities, and create a new bam file
#    ${p2samtools}/samtools view -F 4 -Sb ${out}.ribo.sam > ${out}.ribo.bam & 
#elif [ $stranded == "y" ]
#then
#    ${p2samtools}/samtools view -f 16 ${out}.ribo.sam | awk 'BEGIN {OFS="\n"} {print "@"$1, $10, "+", $11}' > ${out}.nonRibo.fastq
#    ${p2samtools}/samtools view -f 4 ${out}.ribo.sam | awk 'BEGIN {OFS="\n"} {print "@"$1, $10, "+", $11}' > ${out}.nonRibo.fastq &
#    ${p2samtools}/samtools view -f 0 -F 4 -Sb ${out}.ribo.sam > ${out}.ribo.bam &
#fi

#wait

# zip fastq file and delete sam file from first mapping
#gzip ${out}.nonRibo.fastq &
#rm ${out}.ribo.sam

#wait

