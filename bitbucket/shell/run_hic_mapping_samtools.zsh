#!/bin/zsh
# Align Hi-C reads to a genome assembly and assign read pairs to restriction fragments
mini=(-minimap '/nfs/scratch/harmeet.chawla/00_tools_n_scripts/minimap2-2.17_x64-linux/minimap2')
cut=(-cutadapt '/nfs/scratch/harmeet.chawla/anaconda/envs/py35/bin/cutadapt')
sam=(-samtools '/nfs/scratch/harmeet.chawla/00_tools_n_scripts/samtools-1.11/samtools')
bedt=(-bedtools '/nfs/pgsb/commons/apps/bedtools/bedtools2-2.27/bin/bedtools')
novo=(-novosort '/nfs/scratch/harmeet.chawla/00_tools_n_scripts/novocraft/novosort')
bgz=(-bgzip '/nfs/pgsb/commons/apps/bcftools/bcftools-1.3/htslib-1.3/bgzip')
indexsize='50G'
batchmem='5G'

add_dist='
BEGIN{
 OFS="\t"
}

$9 == "-" {
 d1=$2-$12
}

$10 == "-" {
 d2=$5-$15
}

$9 == "+" {
 d1=$13-$3
}

$10 == "+" {
 d2=$16-$6
}

d1 < 0 {
 d1 = 0
}

d2 < 0 {
 d1 = 0
}

$1 == $4 {
 if($12 == $15)
  rd = 0
 else if($13 < $15)
  rd = $15-2-($13+3)
 else
  rd = $12-2-($16+3)
}

$1 != $4 {
 rd=-1
}

{
 l1=d1+$3-$2
 l2=d2+$6-$5
 print $0"\t"d1,d2,l1,l2,l1+l2,rd
}'

dist_stat='
BEGIN{
 OFS="\t"
}

#same fragment
$22 == 0 {
 same++
}

#different seq scaffold
$22 == -1 {
 diffcl++
}

#adj frags
#$22 == 1 {
# adj++
#}

#same seq scaffold
$22 > 1 {
 samecl++
}

END {
 print NR,0+same,0+diffcl,0+samecl
}'

l=(-linker "AAGCTAGCTT") #for HindIII, use GATCGATC for DpnII
m=(-minlen_read 30)
q=(-minq 10)
x=(-maxlen 500)
o=(-onlyfrag 0)

zparseopts -D -K -- -linker:=l -minlen_read:=m -threads:=t \
 -mem:=e -bed:=b -ref:=r -tmp:=p -minq:=q -maxlen:=x -onlyfrag:=o \
 -cutadapt:=cut -minimap:=mini -samtools:=sam \
 -bedtools:=bedt -novosort:=novo -bgzip:=bgz

bgzip=$bgz[2]
novosort=$novo[2]
bedtools=$bedt[2]
samtools=$sam[2]
minimap=$mini[2]
cutadapt=$cut[2]
linker=$l[2]
minlen=$m[2]
threads=$t[2]
mem=$e[2]
ref=$r[2]
tmp=$p[2]
minq=$q[2]
maxlen=$x[2]
bed=$b[2]
onlyfrag=$o[2]
i=$1

if [[ $threads -gt 16 ]]; then
 novosort_threads=16
else
 novosort_threads=$threads
fi

base=$i/$i:t
echo -e "blabla $i"
echo -e "pipip $base"
bam=${base}.bam
minimaperr=${base}_minimap.err
samtoolserr=${base}_samtools.err
sorterr1=${base}_samsort1.err
sorterr2=${base}_samsort2.err
cutadapterr=${base}_cutadapt.err
sorterr0=${base}_samtools0.err

b=$base:t
rgentry="@RG\tID:$b\tPL:ILLUMINA\tPU:$b\tSM:$b"

# Note the find commands. If sample names contain R1 and R2, this will not work.
if [[ $onlyfrag -eq 0 ]]; then
 file=$(find $i | grep _R1 | grep 'f.*q.gz$')
echo -e "poook$i,....$file"
 $cutadapt -f fastq --interleaved -a $linker -A $linker -O 1 -m $minlen 2> $cutadapterr \
   <(find $i | grep _R1 | grep 'f.*q.gz$' | sort | xargs zcat) \
   <(find $i | grep _R2 | grep 'f.*q.gz$' | sort | xargs zcat) \
  | $minimap -ax sr -R $rgentry -t $threads -2 -I $indexsize -K$batchmem $ref /dev/stdin 2> $minimaperr \
  | $samtools fixmate -m - - 2> samtools_fixmate.err \
  | $samtools sort -m 6G -@ 8 2> samtools_sort_pos.err \
  | $samtools markdup -r - - 2> samtools_markdup.err \
  | $samtools sort -@ $novosort_threads -m $mem -n -o $bam 2> samtools_sort_name.err
  
 echo $pipestatus > ${base}_pipestatus.txt
fi

f="${base}_reads_to_fragments.bed.gz"
err="${base}_reads_to_fragments.err"

{
 {
  $samtools view -H $bam
  $samtools view -q 10 -F 3332 $bam \
   | tee >(cut -f 1 | uniq -c | awk '$1 == 2' | wc -l > ${base}_both_mapped_q10.len) \
   | cut -f -16 \
   | awk 'old != $1 {old=$1; printf "\n"$0; next} {printf "\t"$0}' \
   | awk -F '\t' 'NF == 32' \
   | awk '{for(i=1; i<=15; i++) printf $i"\t"; printf $16"\n"$17; for(i=18; i<=32; i++) printf "\t"$i; print ""}'
 } | $samtools view -Su - \
  | $bedtools pairtobed -f 1 -type both -abam - -bedpe -b $bed \
  | awk 'old != $7 {old=$7; printf "\n"$0; next} {printf "\t"$0}' \
  | awk 'NF == 26'  | cut -f 1-13,24-26 | awk "$add_dist" \
  | tee >($bgzip -c > $f) \
        >(wc -l > $f:r.len) \
	>(awk -v maxlen=$maxlen '$21 > maxlen && $1 == $4' \
	   | cut -f 1,2,3,5,6 | awk '{print $5 - $2}' \
	   | awk -v maxlen=$maxlen '$1 <= maxlen' | sort | uniq -c > ${base}_length_dist_PE.txt ) \
        >(awk -v maxlen=$maxlen '$21 > maxlen' | wc -l > ${base}_pe_count.txt) \
	>(awk -v maxlen=$maxlen '$21 <= maxlen' | cut -f 21 | sort | uniq -c > ${base}_length_dist.txt) \
        >(awk -v maxlen=$maxlen '$21 <= maxlen' | awk "$dist_stat" > ${base}_frag_stat.txt) \
 | awk -v maxlen=$maxlen '$21 <= maxlen' \
 | cut -f 11,12,14,15 | awk '$1 != $3 || $2 != $4' \
 | awk '{print $0"\n"$3"\t"$4"\t"$1"\t"$2}' \
 | $bgzip > ${base}_fragment_pairs.tsv.gz
} 2> $err
