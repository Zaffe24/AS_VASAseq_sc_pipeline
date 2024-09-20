#!/bin/bash

if [ $# -ne 4 ]
then
    echo "Please, give (1) input fastq file; (2) output directory; (3) path2trimgalore; (4) path2cutadapt;"
  exit
fi

which cutadapt
which Trim_Galore

file2trim=$1
outdir=$2
path2trimgalore=$3
path2cutadapt=$4
prefix=$(basename "$file2trim")

echo file basename: $prefix
# trim adaptors
${path2trimgalore}/trim_galore --path_to_cutadapt ${path2cutadapt}/cutadapt ${file2trim} -o ${outdir}
sleep 2


# trim homopolymers
${path2cutadapt}/cutadapt -m 15 --trim-n -a "polyG1=GG{5}" -a "polyC1=CC{5}" -a "polyT1=TT{5}" -a "polyA1=AA{5}" -o ${outdir}/${prefix%.fastq.gz}_trimmed_homoATCG.fq.gz  ${outdir}/${prefix%.fastq.gz}_trimmed.fq.gz

sleep 2

file_size=$(stat -c%s ${outdir}/${prefix%.fastq.gz}_trimmed_homoATCG.fq.gz)
min=100
if [ $file_size -gt $min ]; then
    echo -e 'file size in bytes:' "$file_size"
else
    echo -e 'file size in bytes:' "$file_size" 'insufficient'
    exit -1
fi





