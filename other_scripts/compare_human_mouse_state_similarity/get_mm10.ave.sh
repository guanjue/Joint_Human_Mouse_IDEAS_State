#!/bin/bash

#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

source ~/.bashrc

cd /homes1/gxiang/softwares/EpiAlign/Ccode/example/VISION_IDEAS_states
state_bed='S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state.AVE.bed'
while read chr start end st gn
do
echo $st
echo $chr $start $end $st $gn | awk -F ' ' -v OFS='\t' '{print $1,$2,$3,$4,$5}' > 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.bed'
bedtools intersect -a $state_bed -b 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.bed' -wa -u > 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.state.bed'
if [[ $st == "+" ]]
then
cut -f4 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.state.bed' > 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.state.bed.seq'
fi
if [[ $st == "-" ]]
then
cut -f4 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.state.bed' > 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.state.bed.r.seq'
Rscript reverse_seq.R 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.state.bed.r.seq' 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.state.bed.seq'
rm 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.state.bed.r.seq'
fi
Rscript compress_seq.R 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.state.bed.seq' 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.state.bed.Sseq'
rm 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.state.bed'
#rm 'genes/mm10.gene.tss.exp.'$gn'.'$state_bed'.state.bed.seq'
done < genes/mm10.gene.tss.exp.bed

