cd /homes1/gxiang/softwares/EpiAlign/Ccode/example/VISION_IDEAS_states/

### prepare gene set
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.basic.annotation.gff3.gz
#gunzip gencode.v39.basic.annotation.gff3.gz
#cat gencode.v39.basic.annotation.gff3 | awk -F 'gene_type=' '{print $2}' | awk -F ';' '{print $1}' > gene_types.txt
#cat gencode.v39.basic.annotation.gff3 | awk -F 'gene_name=' '{print $2}' | awk -F ';' '{print $1}' > gene_names.txt
#paste gene_types.txt gene_names.txt gencode.v39.basic.annotation.gff3 | awk -F '\t' -v OFS='\t' '{if ($1=="protein_coding" && $5=="gene") print $3,$6,$7,$9,$2}' > hg38.gene.bed
#rm gene_types.txt gene_names.txt

#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gff3.gz
#gunzip gencode.vM25.basic.annotation.gff3.gz
#cat gencode.vM25.basic.annotation.gff3 | awk -F 'gene_type=' '{print $2}' | awk -F ';' '{print $1}' > gene_types.txt
#cat gencode.vM25.basic.annotation.gff3 | awk -F 'gene_name=' '{print $2}' | awk -F ';' '{print $1}' > gene_names.txt
#paste gene_types.txt gene_names.txt gencode.vM25.basic.annotation.gff3 | awk -F '\t' -v OFS='\t' '{if ($1=="protein_coding" && $5=="gene") print $3,$6,$7,$9,$2}' > mm10.gene.bed
#rm gene_types.txt gene_names.txt


### prepare matched ct state bed
#cat S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4, $5, $6,$7, $12, $16,$17, $18, $20, $32,$33, $34, $42, $44, $46  }' > S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state.matched_ct.bed
#cat S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4, $5, $6,$7, $11, $14,$15, $20, $22, $26,$27, $28, $30, $31, $32  }' > S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state.matched_ct.bed


hg38_gene=$1
mm10_gene=$2
hg38_gene_exp_win=$3
mm10_gene_exp_win=$4
hg38_gene_set=$5
mm10_gene_set=$6
hg38_state_set=$7
mm10_state_set=$8


#hg38_gene='GATA1'
#mm10_gene='Gata1'
hg38_gene_set='hg38.gene.bed'
mm10_gene_set='mm10.gene.bed'
hg38_state_set='S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state.matched_ct.bed'
mm10_state_set='S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state.matched_ct.bed'
#hg38_gene_exp_win=50000
#mm10_gene_exp_win=50000

#hg38_gene='HBA1'
#mm10_gene='Hba-a1'
#hg38_gene_set='hg38.gene.bed'
#mm10_gene_set='mm10.gene.bed'
#hg38_state_set='S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state.matched_ct.bed'
#mm10_state_set='S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state.matched_ct.bed'
#hg38_gene_exp_win=100000
#mm10_gene_exp_win=100000

### get query gene
cat $hg38_gene_set | awk -F '\t' -v OFS='\t' -v hg38_gene=$hg38_gene '{if (toupper($5)==toupper(hg38_gene)) print $1,$2,$3,$4,$5}' | awk -F '\t' -v OFS='\t' -v hg38_gene_exp_win=$hg38_gene_exp_win '{if ($4=="+") print $1,$2-hg38_gene_exp_win,$2+hg38_gene_exp_win,$4,$5; else print $1,$3-hg38_gene_exp_win,$3+hg38_gene_exp_win,$4,$5}' | awk -F '\t' -v OFS='\t' '{if ($2<0) print $1,0,$3,$4,$5; else print $0}' > 'hg38.gene.'$hg38_gene'.bed'
cat $mm10_gene_set | awk -F '\t' -v OFS='\t' -v mm10_gene=$mm10_gene '{if (toupper($5)==toupper(mm10_gene)) print $1,$2,$3,$4,$5}' | awk -F '\t' -v OFS='\t' -v mm10_gene_exp_win=$mm10_gene_exp_win '{if ($4=="+") print $1,$2-mm10_gene_exp_win,$2+mm10_gene_exp_win,$4,$5; else print $1,$3-mm10_gene_exp_win,$3+mm10_gene_exp_win,$4,$5}' | awk -F '\t' -v OFS='\t' '{if ($2<0) print $1,0,$3,$4,$5; else print $0}' > 'mm10.gene.'$mm10_gene'.bed'

### intersect state bed
bedtools intersect -a $hg38_state_set -b 'hg38.gene.'$hg38_gene'.bed' -wa -u > 'hg38.gene.'$hg38_gene'.matched_ct.state.bed'
bedtools intersect -a $mm10_state_set -b 'mm10.gene.'$mm10_gene'.bed' -wa -u > 'mm10.gene.'$mm10_gene'.matched_ct.state.bed'

### get state similarity matrix
cat 'hg38.gene.'$hg38_gene'.bed' 'mm10.gene.'$mm10_gene'.bed'
Rscript get_hg38_mm10.cor.heatmap.R 'hg38.gene.'$hg38_gene'.matched_ct.state.bed' 'mm10.gene.'$mm10_gene'.matched_ct.state.bed' $hg38_gene'.cor.heatmap.png' 06a_S3V2_IDEAS_hg38_r3_withHg38Mm10prior.para.modified.para 8 


