# Author: Guanjue Xiang
# Date: 2024-01-01
##############################################################################################################

# time bash /Users/guanjuexiang/Documents/projects/Joint_Human_Mouse_IDEAS_State/other_scripts/cormat_NMF_FDR/Interspec-EpiState-NMF.sh /Users/guanjuexiang/Documents/projects/Joint_Human_Mouse_IDEAS_State/other_scripts/cormat_NMF_FDR/config.info.txt
# Interspec-EpiState-NMF Identifying inter species epigenetic state correlation between human and mouse genes by NMF

##############################################################################################################
# read parameters
config_file=$1
param(){
  cat $config_file | grep -w $1 | cut -f2
}

# global parameters
script_dir=$(param "script_dir")
echo $script_dir
working_dir=$(param "working_dir")
echo $working_dir

# target genes parameters
hg38_gene=$(param "hg38_gene")
mm10_gene=$(param "mm10_gene")
hg38_gene_set=$(param "hg38_gene_set")
mm10_gene_set=$(param "mm10_gene_set")
hg38_gene_exp_win_u=$(param "hg38_gene_exp_win_u")
hg38_gene_exp_win_d=$(param "hg38_gene_exp_win_d")
mm10_gene_exp_win_u=$(param "mm10_gene_exp_win_u")
mm10_gene_exp_win_d=$(param "mm10_gene_exp_win_d")

# epigenetic state parameters
hg38_state_set=$(param "hg38_state_set")
mm10_state_set=$(param "mm10_state_set")
EpigeneticState_meansignal_mat_file=$(param "EpigeneticState_meansignal_mat_file")

# FDR parameters
random_background_gene_num=$(param "random_background_gene_num")
fdr_threshold=$(param "fdr_threshold")
NMF_component_num=$(param "NMF_component_num")

##############################################################################################################

# go to working dir 
cd $working_dir 

# print target gene
echo check input parameters:
echo hg38_gene: $hg38_gene
echo mm10_gene: $mm10_gene
echo hg38_gene_set: $hg38_gene_set
echo mm10_gene_set: $mm10_gene_set
echo hg38_gene_exp_win_u: $hg38_gene_exp_win_u
echo hg38_gene_exp_win_d: $hg38_gene_exp_win_d
echo mm10_gene_exp_win_u: $mm10_gene_exp_win_u
echo mm10_gene_exp_win_d: $mm10_gene_exp_win_d
echo hg38_state_set: $hg38_state_set
echo mm10_state_set: $mm10_state_set
echo EpigeneticState_meansignal_mat_file: $EpigeneticState_meansignal_mat_file
echo random_background_gene_num: $random_background_gene_num
echo fdr_threshold: $fdr_threshold
echo NMF_component_num: $NMF_component_num

##############################################################################################################

# get gene locus coordinate
cat $hg38_gene_set | awk -F '\t' -v OFS='\t' -v hg38_gene=$hg38_gene '{if (toupper($5)==toupper(hg38_gene)) print $1,$2,$3,$4,$5}' | \
awk -F '\t' -v OFS='\t' -v hg38_gene_exp_win_u=$hg38_gene_exp_win_u -v hg38_gene_exp_win_d=$hg38_gene_exp_win_d '{if ($4=="+") print $1,$2-hg38_gene_exp_win_u,$2+hg38_gene_exp_win_d,$4,$5; else print $1,$3-hg38_gene_exp_win_u,$3+hg38_gene_exp_win_d,$4,$5}' | \
awk -F '\t' -v OFS='\t' '{if ($2<0) print $1,0,$3,$4,$5; else print $0}' > 'hg38.gene.'$hg38_gene'.bed'
# mm10 gene
cat $mm10_gene_set | awk -F '\t' -v OFS='\t' -v mm10_gene=$mm10_gene '{if (toupper($5)==toupper(mm10_gene)) print $1,$2,$3,$4,$5}' | \
awk -F '\t' -v OFS='\t' -v mm10_gene_exp_win_u=$mm10_gene_exp_win_u -v mm10_gene_exp_win_d=$mm10_gene_exp_win_d '{if ($4=="+") print $1,$2-mm10_gene_exp_win_u,$2+mm10_gene_exp_win_d,$4,$5; else print $1,$3-mm10_gene_exp_win_u,$3+mm10_gene_exp_win_d,$4,$5}' | \
awk -F '\t' -v OFS='\t' '{if ($2<0) print $1,0,$3,$4,$5; else print $0}' > 'mm10.gene.'$mm10_gene'.bed'

# intersect state bed
bedtools intersect -a $hg38_state_set -b 'hg38.gene.'$hg38_gene'.bed' -wa -u > 'hg38.gene.'$hg38_gene'.matched_ct.state.bed'
bedtools intersect -a $mm10_state_set -b 'mm10.gene.'$mm10_gene'.bed' -wa -u > 'mm10.gene.'$mm10_gene'.matched_ct.state.bed'

# get state similarity matrix
cat 'hg38.gene.'$hg38_gene'.bed' 'mm10.gene.'$mm10_gene'.bed'
Rscript $script_dir/get_hg38_mm10.cor.heatmap.R 'hg38.gene.'$hg38_gene'.matched_ct.state.bed' 'mm10.gene.'$mm10_gene'.matched_ct.state.bed' $hg38_gene'.'$mm10_gene'.cor.heatmap.png' $EpigeneticState_meansignal_mat_file 

##############################################################################################################

# check if $working_dir/randomly_pick_bg_genes_hg38 folder not exist
if [ ! -d $working_dir/randomly_pick_bg_genes_hg38 ]; then
    # select background genes
    # hg38
    Rscript $script_dir/random_select_background_genes.R \
    $working_dir \
    $hg38_gene_set \
    123 \
    $hg38_gene_exp_win_u \
    $hg38_gene_exp_win_d \
    $random_background_gene_num \
    $working_dir/randomly_pick_bg_genes_hg38
    #
    # get epigenetic state correlation matrix for background genes
    # hg38
    ls randomly_pick_bg_genes_hg38/bg.gene.*.bed | awk -F '.' '{print $3}' > randomly_pick_bg_genes_hg38/background.gene.all.txt
    # loop all background genes
    for i in $(cat randomly_pick_bg_genes_hg38/background.gene.all.txt)
    do
    echo $i
    # intersect state bed
    bedtools intersect -a $hg38_state_set -b 'randomly_pick_bg_genes_hg38/bg.gene.'$i'.bed' -wa -u > 'randomly_pick_bg_genes_hg38/bg.gene.'$i'.matched_ct.state.bed'
    # get state similarity matrix
    cat 'mm10.gene.'$mm10_gene'.bed' 'randomly_pick_bg_genes_hg38/bg.gene.'$i'.bed'
    Rscript $script_dir/get_hg38_mm10.cor.heatmap.simple.R 'mm10.gene.'$mm10_gene'.matched_ct.state.bed' 'randomly_pick_bg_genes_hg38/bg.gene.'$i'.matched_ct.state.bed' $mm10_gene'.'$i'.cor.heatmap.png' $EpigeneticState_meansignal_mat_file 
    # rm matched_ct.state.bed
    rm 'randomly_pick_bg_genes_hg38/bg.gene.'$i'.matched_ct.state.bed'
    mv $mm10_gene'.'$i'.cor.heatmap.png.cor.mat.txt' randomly_pick_bg_genes_hg38/
    done
fi

# check if $working_dir/randomly_pick_bg_genes_hg38 folder not exist
if [ ! -d $working_dir/randomly_pick_bg_genes_mm10 ]; then
    # select background genes
    # mm10
    Rscript $script_dir/random_select_background_genes.R \
    $working_dir \
    $mm10_gene_set \
    123 \
    $mm10_gene_exp_win_u \
    $mm10_gene_exp_win_d \
    $random_background_gene_num \
    $working_dir/randomly_pick_bg_genes_mm10
    #
    # get epigenetic state correlation matrix for background genes
    # mm10
    ls randomly_pick_bg_genes_mm10/bg.gene.*.bed | awk -F '.' '{print $3}' > randomly_pick_bg_genes_mm10/background.gene.all.txt
    # loop all background genes
    for i in $(cat randomly_pick_bg_genes_mm10/background.gene.all.txt)
    do
    echo $i
    # intersect state bed
    bedtools intersect -a $mm10_state_set -b 'randomly_pick_bg_genes_mm10/bg.gene.'$i'.bed' -wa -u > 'randomly_pick_bg_genes_mm10/bg.gene.'$i'.matched_ct.state.bed'
    # get state similarity matrix
    cat 'hg38.gene.'$hg38_gene'.bed' 'randomly_pick_bg_genes_mm10/bg.gene.'$i'.bed'
    Rscript $script_dir/get_hg38_mm10.cor.heatmap.simple.R 'hg38.gene.'$hg38_gene'.matched_ct.state.bed' 'randomly_pick_bg_genes_mm10/bg.gene.'$i'.matched_ct.state.bed' $hg38_gene'.'$i'.cor.heatmap.png' $EpigeneticState_meansignal_mat_file 
    # rm matched_ct.state.bed
    rm 'randomly_pick_bg_genes_mm10/bg.gene.'$i'.matched_ct.state.bed'
    mv $hg38_gene'.'$i'.cor.heatmap.png.cor.mat.txt' randomly_pick_bg_genes_mm10/
    done
fi

##############################################################################################################

# get NMF and FDR 
Rscript $script_dir/get_NMF_FDR.R \
$working_dir \
$hg38_gene \
$mm10_gene \
randomly_pick_bg_genes_mm10/background.gene.all.txt \
randomly_pick_bg_genes_hg38/background.gene.all.txt \
$fdr_threshold \
$NMF_component_num

##############################################################################################################

# time bash /Users/guanjuexiang/Documents/projects/Joint_Human_Mouse_IDEAS_State/other_scripts/cormat_NMF_FDR/Interspec-EpiState-NMF.sh /Users/guanjuexiang/Documents/projects/Joint_Human_Mouse_IDEAS_State/other_scripts/cormat_NMF_FDR/input_files/config.info.txt 2> test_run.log.txt
# Platform: MacBook Pro (2021) Apple M1 Processor 16GB RAM
# RunTime: 1941.79s ~32min

