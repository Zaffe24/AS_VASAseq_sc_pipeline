#!/bin/bash

p1=$1
p2=$2
p3=$3
pl=$4
#p1= os.path.join(config["tmp"],config["patient"]),
#p2= config["path_to_plates"],
#p3= config["working_dir"],
#pl=plate,

echo "$@"
echo
echo $(which python)
echo
echo 'plate to be processed:' $pl
echo 
      
mkdir -p $p1/raw_plates/
IFS=';' read -ra array <<< $p2
echo -e "all paths defined:" "${array[@]}"
echo
paths=()

for path in "${array[@]}"; do
    path="${path//$'\n'/}"
    fqs=$(ls ${path}/${pl}_*fastq.gz 2>/dev/null)
    
    if [ "${#fqs}" -gt 2 ]; then
        echo 'found a path:' $path
        paths+=("$path")
    fi
done

echo
echo 'length paths found:' "${#paths[@]}"
echo
echo "${paths[@]}"
echo



if [ "${#paths[@]}" -eq 0 ]; then
    echo -e 'raw paths defined do not contain ' $pl
    exit -1
elif [ "${#paths[@]}" -gt 1 ]; then
    echo -e $pl 'is present in multiple paths'
    exit -1
fi

cd $paths
fqs=$(ls ${pl}_*fastq.gz)

echo
echo 'plate files:' $fqs
echo

if [ "${#fqs[@]}" -eq 0 ]; then
    echo 'files for plate ' $pl 'not found'
    exit -1
fi
#cp -u $fqs -t {params.p1}/raw_plates/
rsync -av --partial --progress $fqs $p1/raw_plates/
cd $p1/raw_plates/

python $p3/scripts/concatenator.py --fqf "$pl" --cbcfile $p3/bc_celseq2.tsv --cbchd 0 --lenumi 6 --umifirst --demux --outdir $p1/raw_fastq/ &
wait
echo
echo demultiplexing for "$pl" concluded
rm $p1/raw_plates/${pl}*