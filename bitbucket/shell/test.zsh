#!/bin/zsh

minimap='/opt/Bio/minimap2/2.17/bin/minimap2'
samtools='/opt/Bio/samtools/1.16.1/bin/samtools'
novosort="/opt/Bio/novocraft/V3.06.05/bin/novosort"
cutadapt='/opt/Bio/cutadapt/3.3/bin/cutadapt'

indexsize='50G'
batchmem='5G'

m=(-minlen_read 30)

zparseopts -D -K -- -minlen_read:=m -threads:=t \
 -mem:=e -ref:=r -tmp:=p  -linker:=l

linker=$l[2]
minlen=$m[2]
threads=$t[2]
mem=$e[2]
ref=$r[2]
tmp=$p[2]
i=$1

if [[ $threads -gt 16 ]]; then
 novosort_threads=16
else
 novosort_threads=$threads
fi

base=$i/$i:t
bam=${base}.bam
minimaperr=${base}_minimap.err
samtoolserr=${base}_samtools.err
sorterr1=${base}_novosort1.err
sorterr2=${base}_novosort2.err
cutadapterr=${base}_cutadapt.err

b=$base:t
rgentry="@RG\tID:$b\tPL:ILLUMINA\tPU:$b\tSM:$b"

$cutadapt --interleaved -a $linker -A $linker \
  -m $minlen 2> $cutadapterr \
  <(find $i | egrep '_R1' | grep 'f.*q$' | sort | xargs cat) \
  <(find $i | egrep '_R2' | grep 'f.*q$' | sort | xargs cat) \
 | $minimap -ax sr -R $rgentry -t $threads -2 -I $indexsize -K$batchmem $ref /dev/stdin 2> $minimaperr \
 | $samtools view -Su /dev/stdin 2> $samtoolserr \
 | $novosort -c $novosort_threads -t $tmp -m $mem --keepTags --md /dev/stdin 2> $sorterr1 \
 | $novosort -c $novosort_threads -t $tmp -m $mem -n -o $bam /dev/stdin 2> $sorterr2

echo $pipestatus > ${base}_pipestatus.txt

$samtools view -q 10 -F 3332 $bam | cut -f -4 \
    | tee >(cut -f 1 | uniq -c | awk '$1 == 2' | wc -l > ${base}_both_mapped_q10.len) \
    | awk 'old != $1 {old=$1; printf "\n"$0; next} {printf "\t"$0}' \
    | awk -F '\t' 'NF == 8' | cut -f 3,4,7,8 \
    | awk '{printf "%s\t%20d\t%s\t%20d\n", $1,1+($2 / 10000)*10000,$3,($4 / 10000)*10000}' \
    | gzip > ${base}_fragment_pairs.tsv.gz 
