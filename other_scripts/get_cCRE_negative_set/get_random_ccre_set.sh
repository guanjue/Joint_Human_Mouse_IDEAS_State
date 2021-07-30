#!/bin/bash

#SBATCH --mem=200G
#SBATCH --cpus-per-task=12

source ~/.bashrc

cd /homes1/gxiang/projects/impute_cistrome/get_neg_set

j=$1
bash -c "RANDOM=$j"

### split orignial cCRE list into 12 peaks
#head -20000 S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.bed > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.1.bed
#head -40000 S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.bed | tail -n+20001 > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.2.bed
#head -60000 S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.bed | tail -n+40001 > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.3.bed
#head -80000 S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.bed | tail -n+60001 > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.4.bed
#head -100000 S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.bed | tail -n+80001 > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.5.bed
#head -120000 S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.bed | tail -n+100001 > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.6.bed
#head -140000 S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.bed | tail -n+120001 > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.7.bed
#head -160000 S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.bed | tail -n+140001 > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.8.bed
#head -180000 S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.bed | tail -n+160001 > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.9.bed
#head -200000 S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.bed | tail -n+180001 > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.10.bed
#head -220000 S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.bed | tail -n+200001 > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.11.bed
#cat S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.bed | tail -n+220001 > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.12.bed


#echo $j
### get random peaks
for i in {1..12}
do
	echo $i
	time bash Negetive_sequence_matched_length_GC.sh S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.$i.bed hg38.chrom.1_22XY.sizes hg38.fa 5 &
done


### pool random peaks
#rm S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.N.bed_matched.$j.bed
#for i in {1..12}
#do
#	cat S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.$i.bed_matched.bed >> S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.N.bed_matched.$j.bed
#done



