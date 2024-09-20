#!/bin/bash

if [ $# -lt 4 ]
       then
     echo ""
     echo "Usage: $0 <human_mouse index folder> <read1.fq.gz> <output_path> <BBMap_path>"
     echo ""
     exit -1
     fi
     

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

#export PATH=/g/korbel/zafferani/bin/:${PATH}

# CMD parameters
HG=${1}
FQ1=${2}
# FQ2=${4}
OP=$(basename "${FQ1}") #EMB-JL-v058_HVMFVBGXN_257.nonRibo.fastq.gz
OUT=${3}   #/scratch/zafferani/P6/cleaned
bbmap=${4}

#pref=$OUT/${OP%.nonRibo*z}
pref=$OUT/${OP%%.nonRibo*.gz}
echo ${OP} ${OUT} ${pref}
# Index the human and mouse reference genomes (required only once)
# bbsplit.sh build=1 ref_x=hs37d5.fa ref_y=mm10.fa

# BBMap
#module load BBMap
"$bbmap"/bbsplit.sh -Xmx32g usejni=t build=1 path=${HG} in=${FQ1} basename="${pref}%_#.fq" minratio=0.5 maxindel=100000 minhits=1 ambiguous2=all local

#echo "${OP%_#.fq}"

# Statistics
cat ${pref}x_1.fq | awk 'NR%4==1' > ${pref}x.reads
cat ${pref}y_1.fq | awk 'NR%4==1' > ${pref}y.reads
#cat AMBIGUOUS_${OP}x_1.fq | awk 'NR%4==1' > AMBIGUOUS_${OP}.reads

MOUSE=`sort ${pref}x.reads ${pref}x.reads ${pref}y.reads | uniq -u | wc -l | cut -f 1`
HUMAN=`sort ${pref}y.reads ${pref}y.reads ${pref}x.reads | uniq -u | wc -l | cut -f 1`
AMBIG=`sort ${pref}x.reads ${pref}y.reads | uniq -d | wc -l | cut -f 1`

#AMBIG=`sort AMBIGUOUS_${OP}.reads | uniq -u | wc -l | cut -f 1`

FRAC=`echo "${MOUSE} / ( ${MOUSE} + ${HUMAN} + ${AMBIG} )" | bc -l`
FRAC2=`echo "${AMBIG} / ( ${MOUSE} + ${HUMAN} + ${AMBIG} )" | bc -l`

if [ -z "$FRAC" ] || [ -z "$FRAC2" ]; then
    echo 'fractions of mouse and/or ambigous reads not calculable'
    exit -1
else
    echo -e 'mouse FRAC:' $FRAC
    echo -e 'ambig FRAC2: ' $FRAC2
fi

rm ${pref}x.reads ${pref}y.reads
rm ${pref}y_1.fq
#rm AMBIGUOUS_${OP}*

mv ${pref}x_1.fq ${pref}.cleaned.fq
sleep 1
gzip ${pref}.cleaned.fq

sleep 5

if [ -s ${pref}.cleaned.fq.gz ]; then
    
    echo "Unique mouse reads" ${MOUSE} > ${pref}.filter.stats
    echo "Unique human reads" ${HUMAN} >> ${pref}.filter.stats
    echo "Ambiguous reads" ${AMBIG} >> ${pref}.filter.stats
    echo "Mouse fraction" ${FRAC} >> ${pref}.filter.stats
    echo "Ambigous fraction" ${FRAC2} >> ${pref}.filter.stats
    echo -e $"${pref##*/}.cleaned \t ${FRAC} \t ${FRAC2}" > ${OUT}/${pref##*/}_fracs.txt

else
    exit -1
fi
