setwd('/Users/guanjuexiang/Documents/projects/analysis/0813_human_mouse_state_compare_heatmap')

mm10_gene_file = 'mm10.gene.bed'
random_seed = 123
tss_up_n200 = 257
tss_down_n200 = 250
bin_size = 200
n_gene = 100
output_folder = 'randomly_pick_bg_genes'
housekeep_genes_list = 'Human_Mouse_Common.csv'

# read in mm10_gene
mm10_gene = read.table(mm10_gene_file, header = F)
#housekeep_genes = read.csv(housekeep_genes_list, header = T, stringsAsFactors = F, sep=';')

# get shared genes
#mm10_gene = mm10_gene_0[mm10_gene_0[,5] %in% housekeep_genes$Mouse,]

set.seed(random_seed)

# randomly pick n genes
mm10_gene_sample = mm10_gene[sample(1:nrow(mm10_gene), n_gene),]
colnames(mm10_gene_sample) = NULL

# get TSS
mm10_gene_sample_tss_pos = mm10_gene_sample[mm10_gene_sample[,4]=='+',c(1,2,2,4,5)]
mm10_gene_sample_tss_neg = mm10_gene_sample[mm10_gene_sample[,4]=='-',c(1,3,3,4,5)]
# pool together
mm10_gene_sample_tss = mm10_gene_sample
mm10_gene_sample_tss[1:nrow(mm10_gene_sample_tss_pos),] = mm10_gene_sample_tss_pos
mm10_gene_sample_tss[(nrow(mm10_gene_sample_tss_pos)+1):nrow(mm10_gene_sample_tss),] = mm10_gene_sample_tss_neg

# expand TSS
mm10_gene_sample_tss_coord = as.matrix(mm10_gene_sample_tss[,2:3])
mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='+',1] = mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='+',1] - tss_up_n200*bin_size
mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='+',2] = mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='+',2] + tss_down_n200*bin_size
mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='-',1] = mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='-',1] - tss_down_n200*bin_size
mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='-',2] = mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='-',2] + tss_up_n200*bin_size
mm10_gene_sample_tss_coord[mm10_gene_sample_tss_coord<0] = 0

# convert to 200bp bin
mm10_gene_sample_tss_coord_bin = mm10_gene_sample_tss_coord/bin_size
mm10_gene_sample_tss_coord_bin = round(mm10_gene_sample_tss_coord_bin)
mm10_gene_sample_tss_coord = mm10_gene_sample_tss_coord_bin*bin_size
mm10_gene_sample_tss_coord[,1] = mm10_gene_sample_tss_coord[,1] + 1

# get TSS
mm10_gene_sample_tss[,2:3] = mm10_gene_sample_tss_coord

# create output folder if not exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# write out TSS
for (i in 1:nrow(mm10_gene_sample_tss)) {
  write.table(mm10_gene_sample_tss[i,], paste0(output_folder, '/', 'mm10.gene.', mm10_gene_sample_tss[i,5], '.bed'), col.names = F, row.names = F, quote = F, sep = '\t')
}


# bash
hg38_chrom_size=/Users/guanjuexiang/Documents/projects/S3V2_IDEAS_ESMP/genomesize/hg38.chrom.1_22XY.sizes
mm10_chrom_size=/Users/guanjuexiang/Documents/projects/S3V2_IDEAS_ESMP/genomesize/mm10.chrom.1_19XY.sizes
script_dir=/Users/guanjuexiang/Documents/projects/Joint_Human_Mouse_IDEAS_State/other_scripts/compare_human_mouse_state_similarity/               
working_dir=/Users/guanjuexiang/Documents/projects/analysis/0813_human_mouse_state_compare_heatmap/

hg38_gene='GATA1'

hg38_gene_set='hg38.gene.bed'
mm10_gene_set='mm10.gene.bed'
hg38_state_set='S3V2_IDEAS_hg38_r3_withHg38Mm10prior.state.matched_ct.bed'
mm10_state_set='S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state.matched_ct.bed'
para_file='06a_S3V2_IDEAS_hg38_r3_withHg38Mm10prior.para.modified.para'

feature_num=8

###
echo $hg38_gene

cd $working_dir
# get all random genes from randomly_pick_bg_genes file names
ls randomly_pick_bg_genes/mm10.gene.*.bed | awk -F '.' '{print $3}' > randomly_pick_bg_genes/mm10.gene.all.txt

# for each gene in mm10_gene, get state similarity matrix
# loop a bash array
for i in $(cat randomly_pick_bg_genes/mm10.gene.all.txt)
do
  echo $i
  # intersect state bed
  bedtools intersect -a $mm10_state_set -b 'randomly_pick_bg_genes/mm10.gene.'$i'.bed' -wa -u > 'randomly_pick_bg_genes/mm10.gene.'$i'.matched_ct.state.bed'
  # get state similarity matrix
  cat 'hg38.gene.'$hg38_gene'.bed' 'mm10.gene.'$i'.bed'
  Rscript $script_dir/get_hg38_mm10.cor.heatmap.simple.R 'hg38.gene.'$hg38_gene'.matched_ct.state.bed' 'randomly_pick_bg_genes/mm10.gene.'$i'.matched_ct.state.bed' $hg38_gene'.'$i'.cor.heatmap.png' $para_file $feature_num 
  # rm matched_ct.state.bed
  rm 'randomly_pick_bg_genes/mm10.gene.'$i'.matched_ct.state.bed'
done



