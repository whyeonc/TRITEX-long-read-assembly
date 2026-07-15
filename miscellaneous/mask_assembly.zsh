#!/bin/zsh

samtools='/opt/Bio/samtools/1.16.1/bin/samtools'
bbduk='/filer-dg/agruppen/DG/mascher/source/BBMap_37.93/bbduk.sh'
bit='/filer-dg/agruppen/DG/mascher/source/faToTwoBit'
tbinfo='/filer-dg/agruppen/DG/mascher/source/twoBitInfo'
bedtools='/opt/Bio/bedtools/2.27.1/bin/bedtools'
kmer='/filer-dg/agruppen/DG/mascher/source/BBMap_37.93/kmercountexact.sh'

s=(-size 31)
m=(-mincount 2)
e=(-mem '100G')

zparseopts -D -K -- -size:=s -fasta:=f -mincount:=m -out:=d -mem:=e

fa=$f[2]
size=$s[2]
mincount=$m[2]
out=$d[2]
mem=$e[2]

mkdir -p $out
base=$out/${fa:t:r}

$kmer -Xmx$mem mincount=$mincount in=$fa \
 out=${base}_kmercount.txt k=$size \
  > ${base}_kmercount.out 2> ${base}_kmercount.err 

$bbduk kmask='N' k=$size -Xmx$mem in=$fa ref=${base}_kmercount.txt \
 out=${base}_masked.fasta > ${base}_bbduk.out 2> ${base}_bbduk.err

$samtools faidx ${base}_masked.fasta 

fa="${base}_masked.fasta"
$bit -long $fa ${fa:r}.bit
$tbinfo -nBed ${fa:r}.bit ${fa:r}_gaps.bed

$bedtools complement -g =(cut -f -2 $fa.fai) -i ${fa:r}_gaps.bed > ${fa:r}_noGaps.bed
