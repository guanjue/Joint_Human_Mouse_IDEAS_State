Working_folder='/Users/guanjuexiang/Documents/projects/analysis/04_05_2022_coe/coe_analysis'
state_coe_score_Script_folder='/Users/guanjuexiang/Documents/projects/Joint_Human_House_IDEAS_State/other_scripts/state_coe_score/'



############################################################################################################
### Start: go to working folder
############################################################################################################
cd $Working_folder
############################################################################################################
### download data
############################################################################################################
wget https://usevision.org/data/hg38/RNA/Oct2021/tpmFeb21_v3.tab
wget https://usevision.org/data/hg38/RNA/old/HumanVISION_RNAseq_hg38_genes_tpm.txt

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz
gunzip gencode.v22.annotation.gtf.gz

wget https://usevision.org/data/hg38/IDEASstates/ideasJointMay2021/S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.bed
wget https://usevision.org/data/hg38/IDEASstates/ideasJointMay2021/S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state
############################################################################################################



############################################################################################################
### preprocess downloaded data
############################################################################################################
### get RNA-seq bed file & matrix
mv HumanVISION_RNAseq_hg38_genes_tpm.txt HumanVISION_RNAseq_hg38_genes_tpm.old.txt
time Rscript $state_coe_score_Script_folder'/coe_human/get_new_RNA_mat.R'
cut -f1,2,3,4 HumanVISION_RNAseq_hg38_genes_tpm.txt | tail -n+2 | sort -k1,1 -k2,2n > RNA_gene.bed
cat HumanVISION_RNAseq_hg38_genes_tpm.txt | head -1 > HumanVISION_RNAseq_hg38_genes_tpm.idsort.header.txt
tail -n+2 HumanVISION_RNAseq_hg38_genes_tpm.txt | sort -k4,4 > HumanVISION_RNAseq_hg38_genes_tpm.idsort.txt

### get gene annotations
cut -f7,9,10 gencode.v22.annotation.gtf > annotation1.txt
cat annotation1.txt | awk -F ' ' -v OFS='\t' '{if ($4=="gene_type") print $1,$3,$5}' > annotation2.txt
cat annotation2.txt | awk -F '"' -v OFS='\t' '{print $2,$4,$1}' | sort -k1,1 > annotation3.idsort.od.txt
time Rscript $state_coe_score_Script_folder'/coe_human/get_new_gene_annotation.R'

### get protein coding genes TPM mat
paste annotation3.idsort.txt HumanVISION_RNAseq_hg38_genes_tpm.idsort.txt \
| awk -F '\t' -v OFS='\t' '{if ($2=="protein_coding") print $5,$6,$7,$2,$3,$8, $9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47,$48,$49,$50,$51,$52}' \
> HumanVISION_RNAseq_hg38_genes_tpm.idsort.protein_coding.txt
### get protein coding genes bed
cut -f1,2,3,5 HumanVISION_RNAseq_hg38_genes_tpm.idsort.protein_coding.txt > HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.bed

### get All genes Proximal regions (+/- 1Kb) bed
paste annotation3.idsort.txt HumanVISION_RNAseq_hg38_genes_tpm.idsort.txt \
| awk -F '\t' -v OFS='\t' '{print $5,$6,$7,$2,$3,$8}' | cut -f1,2,3,5\
| awk -F '\t' -v OFS='\t' -v exp_win=1000 '{if (($2-exp_win>0) && ($4=="+")) print $1,$2-exp_win, $2+exp_win; if (($3-exp_win>0) && ($4=="-")) print $1,$3-exp_win, $3+exp_win; if (($2-exp_win<0) && ($4=="+")) print $1,0, $2+exp_win; if (($3-exp_win<0) && ($4=="-")) print $1,0, $3+exp_win}' \
> HumanVISION_RNAseq_hg38_gene.idsort.all.Nkbupdownexp.bed

### get Protein Coding gene Proximal regions (+/- 1Kb) bed
cat HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.bed \
| awk -F '\t' -v OFS='\t' -v exp_win=1000 '{if (($2-exp_win>0) && ($4=="+")) print $1,$2-exp_win, $2+exp_win; if (($3-exp_win>0) && ($4=="-")) print $1,$3-exp_win, $3+exp_win; if (($2-exp_win<0) && ($4=="+")) print $1,0, $2+exp_win; if (($3-exp_win<0) && ($4=="-")) print $1,0, $3+exp_win}' \
> HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.bed

### get Protein Coding gene Distal regions (+/- 50Kb) bed
cat HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.bed \
| awk -F '\t' -v OFS='\t' -v exp_win=50000 '{if (($2-exp_win>0) && ($4=="+")) print $1,$2-exp_win, $2+exp_win; if (($3-exp_win>0) && ($4=="-")) print $1,$3-exp_win, $3+exp_win; if (($2-exp_win<0) && ($4=="+")) print $1,0, $2+exp_win; if (($3-exp_win<0) && ($4=="-")) print $1,0, $3+exp_win}' \
> HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.bed

### add cCRE ids
bedtools intersect -a S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.bed -b HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.bed -v > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.noProximal.bed
time Rscript $state_coe_score_Script_folder'/coe_human/add_id.R' S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.bed S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.withid.bed

### get cCREs intersect with All genes' Proximal regions
bedtools intersect -a S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.withid.bed -b HumanVISION_RNAseq_hg38_gene.idsort.all.Nkbupdownexp.bed -c > S3V2_IDEAS_hg38_ccre2.cCRE.M.withid.atProximal.bed
############################################################################################################



############################################################################################################
### Get State coverage matrices for genes
############################################################################################################
### convert IDEAS outputs .state file to .bed file
tail -n+2 S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state \
| awk -F ' ' -v OFS='\t' '{print $2,$3,$4, $0}' > J_IDEAS.state.bed

### only save the state bin intersect with cCRE
bedtools intersect -a J_IDEAS.state.bed -b S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.bed -wa -u \
> J_IDEAS.state.ccre.bed
cut -f4 J_IDEAS.state.ccre.bed > J_IDEAS.state.inccre

### only save the state bin intersect with Protein Coding genes Proximal regions
bedtools intersect -a J_IDEAS.state.bed -b HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.bed -wa -u \
> J_IDEAS.state.gene.bed
cut -f4 J_IDEAS.state.gene.bed > J_IDEAS.state.ingene

### get Proximal regions' states coverages matrices (25 states' gene-by-celltype matrices)
for i in {0..24}
do
   echo $i
   cp HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.bed HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.S$i.bed
   for j in {5..47}
   do
      echo $j
      cat J_IDEAS.state.ingene \
      | awk -F ' ' -v OFS='\t' -v c=$j -v s=$i '{if ($c==s) print $2,$3,$4}' > ct$j.s$i.bed
      echo 'NO filter cCRE'
      bedtools intersect -a HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.bed -b ct$j.s$i.bed -c > ct$j.s$i.bed.count.bed
      cut -f4 ct$j.s$i.bed.count.bed > ct$j.s$i.bed.count.txt
      paste HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.S$i.bed ct$j.s$i.bed.count.txt > HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.S$i.bed.tmp \
      && mv HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.S$i.bed.tmp HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.S$i.bed
      rm ct$j.s$i.bed.count.txt ct$j.s$i.bed.count.bed ct$j.s$i.bed
   done
done

### get Distal regions' states coverages matrices (25 states' gene-by-celltype matrices)
for i in {0..24}
do
   echo $i
   cp HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.bed HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.S$i.bed
   for j in {5..47}
   do
      echo $j
      cat J_IDEAS.state.inccre \
      | awk -F ' ' -v OFS='\t' -v c=$j -v s=$i '{if ($c==s) print $2,$3,$4}' > ct$j.s$i.ccre.bed
      echo 'filter cCRE'
      bedtools intersect -a HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.bed -b ct$j.s$i.ccre.bed -c > ct$j.s$i.bed.count.bed
      cut -f4 ct$j.s$i.bed.count.bed > ct$j.s$i.bed.count.txt
      paste HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.S$i.bed ct$j.s$i.bed.count.txt > HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.S$i.bed.tmp \
      && mv HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.S$i.bed.tmp HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.S$i.bed
      rm ct$j.s$i.bed.count.txt ct$j.s$i.bed.count.bed ct$j.s$i.ccre.bed
   done
done

### add id to Protein Coding gene's Distal regions
time Rscript $state_coe_score_Script_folder'/coe_human/add_id.R' HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.bed HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.withid.bed 
### sort by gene-id
sort -k1,1 -k2,2n HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.withid.bed > HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.withid.sort.bed
### get the cCREs' assignments for all Protein Coding genes
bedtools map -a HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.withid.sort.bed -b S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.withid.bed -c 4 -o distinct > HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.withccreid.bed
### sort by gene-id
sort -k4,4n HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.withccreid.bed > HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.withccreid.bed.tmp
### rename
mv HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.withccreid.bed.tmp HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.withccreid.bed
############################################################################################################



############################################################################################################
### Get State coverage matrices for cCREs
############################################################################################################
### get cCREs states coverages matrices (25 states' cCRE-by-celltype matrices)
for i in {0..24}
do
   echo $i
   cp S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.withid.bed S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.withid.S$i.mat.txt
   for j in {5..47}
   do
      echo $j
      cat J_IDEAS.state.inccre \
      | awk -F ' ' -v OFS='\t' -v c=$j -v s=$i '{if ($c==s) print $2,$3,$4}' > ct$j.s$i.bed
      ### NO filter cCRE
      bedtools intersect -a S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.withid.bed -b ct$j.s$i.bed -c > ct$j.s$i.bed.count.bed
      cut -f5 ct$j.s$i.bed.count.bed > ct$j.s$i.bed.count.txt
      paste S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.withid.S$i.mat.txt ct$j.s$i.bed.count.txt > S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.withid.S$i.mat.txt.tmp \
      && mv S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.withid.S$i.mat.txt.tmp S3V2_IDEAS_hg38_ccre2.cCRE.M.notall0.rmallNEU.withid.S$i.mat.txt
      rm ct$j.s$i.bed.count.txt ct$j.s$i.bed.count.bed ct$j.s$i.bed
   done
done
############################################################################################################



############################################################################################################
### Get Beta coefficients using PCA-lm model with leave-one-out
############################################################################################################
declare -a common_ct=("B" "CLP" "CMP" "EOS" "ERY" "GMP" "LSK" "MEP" "MK" "MONc" "MONp" "MPP" "NEU" "NK" "CD4" "CD8")
## loop through the common_ct list
for ct in "${common_ct[@]}"
do
   echo "$ct"
   time Rscript $state_coe_score_Script_folder'/coe_human/get_state_Beta_coefficients.human.R' $ct
done

### get Beta coefficients for Human by taking ave beta cross all leave-one-out runs
Rscript $state_coe_score_Script_folder'/coe_human/get_ave.R'
############################################################################################################



############################################################################################################
### Get average Beta coefficients between Human and Mouse (Needs to be run after Mouse analysis is Completed)
############################################################################################################
### get mouse human ave beta
Rscript $state_coe_score_Script_folder'/coe_human/get_ave_cross_species.R'
############################################################################################################



############################################################################################################
### Get cCRE esRP mat in both species
############################################################################################################
### get cCRE esRP mat in both species
Rscript $state_coe_score_Script_folder'/get_cCRE_coe_mat.R'



