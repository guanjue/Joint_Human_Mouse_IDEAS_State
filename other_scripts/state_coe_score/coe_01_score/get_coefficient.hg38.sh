cd /gpfs/group/yzz2/default/legacy/group/projects/vision_joint_human_mouse/coefficients_human

wget http://usevision.org/data/hg38/HumanVISION_RNAseq_hg38_genes_tpm.txt
cut -f1,2,3,4 HumanVISION_RNAseq_hg38_genes_tpm.txt | tail -n+2 | sort -k1,1 -k2,2n > RNA_gene.bed
tail -n+2 HumanVISION_RNAseq_hg38_genes_tpm.txt | sort -k4,4 > HumanVISION_RNAseq_hg38_genes_tpm.idsort.txt

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz
gunzip gencode.v22.annotation.gtf.gz
cut -f7,9,10 gencode.v22.annotation.gtf > annotation1.txt
cat annotation1.txt | awk -F ' ' -v OFS='\t' '{if ($4=="gene_type") print $1,$3,$5}' > annotation2.txt
cat annotation2.txt | awk -F '"' -v OFS='\t' '{print $2,$4,$1}' | sort -k1,1 > annotation3.idsort.txt


paste annotation3.idsort.txt HumanVISION_RNAseq_hg38_genes_tpm.idsort.txt | awk -F '\t' -v OFS='\t' '{if ($2=="protein_coding") print $5,$6,$7,$2,$3,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39}' > HumanVISION_RNAseq_hg38_genes_tpm.idsort.protein_coding.txt

cut -f1,2,3,5 HumanVISION_RNAseq_hg38_genes_tpm.idsort.protein_coding.txt > HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.bed

cat HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.bed \
| awk -F '\t' -v OFS='\t' -v exp_win=2500 '{if (($2-exp_win>0) && ($4=="+")) print $1,$2-exp_win, $2+exp_win; if (($3-exp_win>0) && ($4=="-")) print $1,$3-exp_win, $3+exp_win; if (($2-exp_win<0) && ($4=="+")) print $1,0, $2+exp_win; if (($3-exp_win<0) && ($4=="-")) print $1,0, $3+exp_win}' \
> HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.bed

cat HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.bed \
| awk -F '\t' -v OFS='\t' -v exp_win=50000 '{if (($2-exp_win>0) && ($4=="+")) print $1,$2-exp_win, $2+exp_win; if (($3-exp_win>0) && ($4=="-")) print $1,$3-exp_win, $3+exp_win; if (($2-exp_win<0) && ($4=="+")) print $1,0, $2+exp_win; if (($3-exp_win<0) && ($4=="-")) print $1,0, $3+exp_win}' \
> HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.bed



for i in {0..23}
do
echo $i
cp HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.bed HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.S$i.bed
for j in {5..43}
do
echo $j
tail -n+2 ~/scratch/S3V2norm_compare/human_vision_S3V2/hg38bp0402_forJoint_r60_run50_MP_IDEAS_output_A1_rerun/hg38bp0402_forJoint_r60_run50_MP.state \
| awk -F ' ' -v OFS='\t' -v c=$j -v s=$i '{if ($c==s) print $2,$3,$4}' > ct$j.s$i.bed
bedtools intersect -a ct$j.s$i.bed -b /storage/home/gzx103/group/projects/vision_human/bp24_v2_cCRE_max1_IDEAS_output_A1/bp24_v2_cCRE_max1.cCRE.M.bed -wa -u > ct$j.s$i.ccre.bed
bedtools intersect -a HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.bed -b ct$j.s$i.ccre.bed -c > ct$j.s$i.bed.count.bed
cut -f4 ct$j.s$i.bed.count.bed > ct$j.s$i.bed.count.txt
paste HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.S$i.bed ct$j.s$i.bed.count.txt > HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.S$i.bed.tmp \
&& mv HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.S$i.bed.tmp HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.S$i.bed
rm ct$j.s$i.bed.count.txt ct$j.s$i.bed.count.bed ct$j.s$i.bed
done
done

for i in {0..23}
do
echo $i
cp HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.bed HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.S$i.bed
for j in {5..43}
do
echo $j
tail -n+2 ~/scratch/S3V2norm_compare/human_vision_S3V2/hg38bp0402_forJoint_r60_run50_MP_IDEAS_output_A1_rerun/hg38bp0402_forJoint_r60_run50_MP.state \
| awk -F ' ' -v OFS='\t' -v c=$j -v s=$i '{if ($c==s) print $2,$3,$4}' > ct$j.s$i.bed
bedtools intersect -a ct$j.s$i.bed -b /storage/home/gzx103/group/projects/vision_human/bp24_v2_cCRE_max1_IDEAS_output_A1/bp24_v2_cCRE_max1.cCRE.M.bed -wa -u > ct$j.s$i.ccre.bed
bedtools intersect -a HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.bed -b ct$j.s$i.ccre.bed -c > ct$j.s$i.bed.count.bed
cut -f4 ct$j.s$i.bed.count.bed > ct$j.s$i.bed.count.txt
paste HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.S$i.bed ct$j.s$i.bed.count.txt > HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.S$i.bed.tmp \
&& mv HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.S$i.bed.tmp HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.S$i.bed
rm ct$j.s$i.bed.count.txt ct$j.s$i.bed.count.bed ct$j.s$i.bed
done
done



time Rscript get_eRP_pcareg.R







