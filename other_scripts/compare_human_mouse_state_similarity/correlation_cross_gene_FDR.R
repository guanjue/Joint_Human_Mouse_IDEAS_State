setwd('/Users/guanjuexiang/Documents/projects/analysis/0813_human_mouse_state_compare_heatmap')

library(LSD)
smooth_cor_mat = function(cor_mat, spar = 0.1) {
  cor_mat_smooth_col = t(apply(cor_mat, 1, function(x) smooth.spline(1:ncol(cor_mat), x, spar = spar)$y))
  cor_mat_smooth_row = apply(cor_mat, 2, function(x) smooth.spline(1:nrow(cor_mat), x, spar = spar)$y)
  cor_mat_smooth = (cor_mat_smooth_col + cor_mat_smooth_row)/2
  return(cor_mat_smooth)
}

plot_heatscatter = function(cor_mat, cor_mat_local, outputname){
    # plot heatscatter of cor_mat
    cor_mat_heatscatter_x = c()
    cor_mat_heatscatter_y = c()
    for (i in 1:nrow(cor_mat)) {
        cor_mat_heatscatter_x = c(cor_mat_heatscatter_x, 1:ncol(cor_mat))
        cor_mat_heatscatter_y = c(cor_mat_heatscatter_y, as.numeric(cor_mat[i,]))
    }
    # plot cor_mat along the columns
    png(outputname, width=1000, height=400)
    heatscatter(cor_mat_heatscatter_x, cor_mat_heatscatter_y, pch = 20, cex = 1, xlab = 'Column', ylab = 'Correlation', ylim=c(min(0),max((1))))
    for (i in 1:nrow(cor_mat_local)) {
        points((1:ncol(cor_mat_local))[cor_mat_local[i,]!=0], (cor_mat_local[i,][cor_mat_local[i,]!=0]), cex = 1, col='black')
    }
    dev.off()
}

hg38_gene = 'GATA1'
mm10_gene = 'Gata1'

mm10_gene_bg0 = 'Cd4'
mm10_gene_bg0 = 'Slc4a1'
mm10_gene_bg0 = 'Rps19'

mm10_gene_bg0 = 'Rps19'


input_correlation_mat = paste0('GATA1.Gata1.cor.heatmap.png.cor.mat.RAW_data_for_plot.txt')
input_correlation_mat_bg = paste0(hg38_gene, '.', mm10_gene_bg0, '.cor.heatmap.png.cor.mat.txt')

fdr_thresh = 1e-2
random_picked_bg_genes_file = 'randomly_pick_bg_genes/mm10.gene.all.txt'
smooth_spar = 0.01

# read in correlation matrix
cor_mat = read.table(input_correlation_mat, header = F)
# read in correlation matrix
cor_mat_bg = read.table(input_correlation_mat_bg, header = F)

# smooth cor_mat matrix along the columns for each row
cor_mat_smooth = smooth_cor_mat(cor_mat, spar = smooth_spar)
cor_mat_bg_smooth = smooth_cor_mat(cor_mat_bg, spar = smooth_spar)


#############################################################################################################
# RAW within-gene adjustment
# foreground
get_2r_fdr = function(cor_mat_smooth, fdr_thresh){
    cor_mat_smooth_zcol_raw = apply(cor_mat_smooth, 2, function(x) (x-mean(x))/sd(x))
    cor_mat_zpfdr_raw = apply(cor_mat_smooth_zcol_raw, 2, function(x) (pnorm(x, lower.tail=F)))
    cor_mat_zpfdr_raw_mean = c()
    cor_mat_zpfdr_raw_sd = c()
    p_thresh = 1e-2
    for (j in 1:dim(cor_mat_smooth)[2]) {
        if (sum(cor_mat_zpfdr_raw[,j]>=fdr_thresh)>1){
            cor_mat_zpfdr_raw_mean = c(cor_mat_zpfdr_raw_mean, mean(cor_mat_smooth[cor_mat_zpfdr_raw[,j]>=p_thresh,j]))
            cor_mat_zpfdr_raw_sd = c(cor_mat_zpfdr_raw_sd, sd(cor_mat_smooth[cor_mat_zpfdr_raw[,j]>=p_thresh,j]))
        } else {
            cor_mat_zpfdr_raw_mean = c(cor_mat_zpfdr_raw_mean, mean(cor_mat_smooth[,j]))
            cor_mat_zpfdr_raw_sd = c(cor_mat_zpfdr_raw_sd, sd(cor_mat_smooth[,j]))
        }
    }
    cor_mat_zpfdr_raw_2 = cor_mat_smooth
    for (j in 1:dim(cor_mat_smooth)[2]) {
        cor_mat_zpfdr_raw_2[,j] = p.adjust(pnorm((cor_mat_smooth[,j] - cor_mat_zpfdr_raw_mean[j]) / cor_mat_zpfdr_raw_sd[j], lower.tail=F), 'fdr')
    }
    cor_mat_local_raw = cor_mat_smooth * (cor_mat_zpfdr_raw_2<fdr_thresh)
    return(cor_mat_local_raw)
}

# foreground
cor_mat_local_raw = get_2r_fdr(cor_mat_smooth, fdr_thresh)
print(sum(cor_mat_local_raw!=0))
print(sum(apply(cor_mat_local_raw, 2, max)>0))
# write cor_mat_local
write.table(cor_mat_local_raw, paste0(hg38_gene, '.', mm10_gene, '.cor.heatmap.png.cor.mat.local.raw.txt'), row.names = F, col.names = F, quote = F, sep = '\t')

# background
cor_mat_bg_local_raw = get_2r_fdr(cor_mat_bg_smooth, fdr_thresh)
print(sum(cor_mat_bg_local_raw!=0))
print(sum(apply(cor_mat_bg_local_raw, 2, max)>0))
# write cor_mat_local
write.table(cor_mat_bg_local_raw, paste0(hg38_gene, '.', mm10_gene_bg0, '.cor.heatmap.png.cor.mat.local.bg.raw.txt'), row.names = F, col.names = F, quote = F, sep = '\t')

print(sum(cor_mat_bg_local_raw!=0)/(sum(cor_mat_bg_local_raw!=0)+sum(cor_mat_local_raw!=0)))
c('Cd4',0.4943951)
c('Rps19',0.2548424)
c('Slc4a1',0.4008246)


# plot correlation manhattan plot
plot_heatscatter(cor_mat_smooth, cor_mat_local_raw, paste0(hg38_gene, '.', mm10_gene, '.cor_column.cor.r.raw.png'))
plot_heatscatter(cor_mat_bg_smooth, cor_mat_bg_local_raw, paste0(hg38_gene, '.', mm10_gene_bg0, '.cor_column.cor.r.bg.raw.png'))

#############################################################################################################



# read random_pick_bg_genes
random_picked_bg_genes = read.table(random_picked_bg_genes_file, header = F)

# read cor_mat for random picked bg genes
cor_mat_bg_r_smooth_mat = matrix(-1000, nrow=dim(random_picked_bg_genes)[1]*dim(cor_mat_bg)[1], ncol=dim(cor_mat_bg)[2])
cor_mat_bg_r_smooth_mat_all = matrix(-1000, nrow=dim(random_picked_bg_genes)[1]*dim(cor_mat_bg)[1], ncol=dim(cor_mat_bg)[2])
cor_mat_bg_r_smooth_mat_start = rep(1, dim(cor_mat_bg)[2])
#random_picked_bg_genes = rbind(c('Cd4', 'Cd4'))

for (i in 1:nrow(random_picked_bg_genes)) {
#for (i in 1:2) {
  mm10_gene_bg = random_picked_bg_genes[i,1]
  print(mm10_gene_bg)
  #mm10_gene_bg = random_picked_bg_genes[1,1]
  #mm10_gene_bg = 'Gata1'
  input_correlation_mat_bg = paste0(hg38_gene, '.', mm10_gene_bg, '.cor.heatmap.png.cor.mat.txt')
  if (file.exists(input_correlation_mat_bg)) {
  cor_mat_bg_r = read.table(input_correlation_mat_bg, header = F)
  cor_mat_bg_r_smooth = smooth_cor_mat(cor_mat_bg_r, spar = smooth_spar)
  # add to cor_mat_bg_r_smooth_mat_all
  cor_mat_bg_r_smooth_mat_all[(1:dim(cor_mat_bg_r_smooth)[1])+(i-1)*dim(cor_mat_bg_r_smooth)[1],] = cor_mat_bg_r_smooth
  # get fdr top genes
  for (j in 1:dim(cor_mat_bg_r_smooth)[2]){
    cor_mat_bg_r_smooth_p_j = pnorm((cor_mat_bg_r_smooth[,j]-mean(cor_mat_bg_r_smooth[,j]))/sd(cor_mat_bg_r_smooth[,j]), lower.tail=F)
    cor_mat_bg_r_smooth_fdr_j = p.adjust(cor_mat_bg_r_smooth_p_j, method = 'fdr')
    # get cor_mat high signals based zpfdr
    #print(sum(cor_mat_bg_r_smooth_fdr_j<fdr_thresh))
    if (sum(cor_mat_bg_r_smooth_fdr_j<fdr_thresh)>1){
      cor_mat_bg_r_smooth_mat[cor_mat_bg_r_smooth_mat_start[j]:(cor_mat_bg_r_smooth_mat_start[j]+sum(cor_mat_bg_r_smooth_fdr_j<fdr_thresh)-1),j] = cor_mat_bg_r_smooth[cor_mat_bg_r_smooth_fdr_j<fdr_thresh,j]
      cor_mat_bg_r_smooth_mat_start[j] = cor_mat_bg_r_smooth_mat_start[j] + sum(cor_mat_bg_r_smooth_fdr_j<fdr_thresh) 
    #} else if (sum(cor_mat_bg_r_smooth_p_j<p_thresh)>1){
    #    cor_mat_bg_r_smooth_mat[cor_mat_bg_r_smooth_mat_start[j]:(cor_mat_bg_r_smooth_mat_start[j]+sum(cor_mat_bg_r_smooth_p_j<p_thresh)-1),j] = cor_mat_bg_r_smooth[cor_mat_bg_r_smooth_p_j<p_thresh,j]
    #    cor_mat_bg_r_smooth_mat_start[j] = cor_mat_bg_r_smooth_mat_start[j] + sum(cor_mat_bg_r_smooth_p_j<p_thresh)
    } else{
        cor_mat_bg_r_smooth_mat[cor_mat_bg_r_smooth_mat_start[j]:(cor_mat_bg_r_smooth_mat_start[j]+sum(cor_mat_bg_r_smooth_fdr_j<=1)-1),j] = cor_mat_bg_r_smooth[cor_mat_bg_r_smooth_fdr_j<=1,j]
        cor_mat_bg_r_smooth_mat_start[j] = cor_mat_bg_r_smooth_mat_start[j] + sum(cor_mat_bg_r_smooth_fdr_j<=1)
    }
  }
  print(sum(cor_mat_bg_r_smooth_mat_start!=1))
}
}

# get mean and sd for each column based on random picked false discovery bins
cor_mat_bg_r_smooth_mean = apply(cor_mat_bg_r_smooth_mat, 2, function(x) mean(x[x!=-1000]))
cor_mat_bg_r_smooth_sd = apply(cor_mat_bg_r_smooth_mat, 2, function(x) sd(x[x!=-1000]))
cor_mat_bg_r_smooth_mean_all = apply(cor_mat_bg_r_smooth_mat_all, 2, function(x) mean(x[x!=-1000]))
cor_mat_bg_r_smooth_sd_all = apply(cor_mat_bg_r_smooth_mat_all, 2, function(x) sd(x[x!=-1000]))
cor_mat_bg_r_smooth_mean[is.na(cor_mat_bg_r_smooth_mean)] = cor_mat_bg_r_smooth_mean_all[is.na(cor_mat_bg_r_smooth_mean)]
cor_mat_bg_r_smooth_sd[is.na(cor_mat_bg_r_smooth_sd)] = cor_mat_bg_r_smooth_sd_all[is.na(cor_mat_bg_r_smooth_sd)]

# FDR along the columns
cor_mat_zcol = matrix(0, nrow=dim(cor_mat)[1], ncol=dim(cor_mat)[2])
non0_col = c()
for (j in 1:dim(cor_mat)[2]) {
  # get cor_mat high signals based zpfdr
  if (cor_mat_bg_r_smooth_mean[j]!=0){
    cor_mat_zcol[,j] = (cor_mat_smooth[,j] - cor_mat_bg_r_smooth_mean[j]) / cor_mat_bg_r_smooth_sd[j]
  } else {
    cor_mat_zcol[,j] = (cor_mat_smooth[,j] - mean(cor_mat_smooth[,j])) / sd(cor_mat_smooth[,j])
  }
}

cor_mat_zpfdr = apply(cor_mat_zcol, 2, function(x) p.adjust(pnorm(x, lower.tail=F), method = 'fdr'))
cor_mat_local = cor_mat_smooth * (cor_mat_zpfdr<fdr_thresh)
print(sum(cor_mat_zpfdr<fdr_thresh))
print(sum(apply(cor_mat_local, 2, max)>0))
# write cor_mat_local
write.table(cor_mat_local, paste0(hg38_gene, '.', mm10_gene, '.cor.heatmap.png.cor.mat.local.txt'), row.names = F, col.names = F, quote = F, sep = '\t')

# negative control Cd4 gene
cor_mat_bg_zcol = matrix(0, nrow=dim(cor_mat)[1], ncol=dim(cor_mat)[2])
non0_col = c()
for (j in 1:dim(cor_mat)[2]) {
  # get cor_mat high signals based zpfdr
  if (cor_mat_bg_r_smooth_mean[j]!=0){
    cor_mat_bg_zcol[,j] = (cor_mat_bg_smooth[,j] - cor_mat_bg_r_smooth_mean[j]) / cor_mat_bg_r_smooth_sd[j]
  }else {
    cor_mat_bg_zcol[,j] = (cor_mat_bg_smooth[,j] - mean(cor_mat_bg_smooth[,j])) / sd(cor_mat_bg_smooth[,j])
  }
}

cor_mat_bg_zpfdr = apply(cor_mat_bg_zcol, 2, function(x) p.adjust(pnorm(x, lower.tail=F), method = 'fdr'))
cor_mat_bg_local = cor_mat_bg_smooth * (cor_mat_bg_zpfdr<fdr_thresh)
print(sum(cor_mat_bg_zpfdr<fdr_thresh))
print(sum(apply(cor_mat_bg_local, 2, max)>0))
# write cor_mat_local
write.table(cor_mat_bg_local, paste0(hg38_gene, '.', mm10_gene_bg0, '.cor.heatmap.png.cor.mat.local.bg.txt'), row.names = F, col.names = F, quote = F, sep = '\t')
print(sum(cor_mat_bg_zpfdr<fdr_thresh)/(sum(cor_mat_bg_zpfdr<fdr_thresh)+sum(cor_mat_zpfdr<fdr_thresh)))
c('Cd4',0.04179168)
c('Rps19',0.2713494)
c('Slc4a1',0.375733)


# plot correlation manhattan plot
plot_heatscatter(cor_mat_smooth, cor_mat_local, paste0(hg38_gene, '.', mm10_gene, '.cor_column.cor.r.png'))
plot_heatscatter(cor_mat_bg_smooth, cor_mat_bg_local, paste0(hg38_gene, '.', mm10_gene_bg0, '.cor_column.cor.r.bg.png'))


#############################################################################################################

FP_n = c()
FP_n_naive = c()
all_bgn = c()

for (i in 1:nrow(random_picked_bg_genes)) {
    mm10_gene_bg_i = random_picked_bg_genes[i,1]
# get bg gene i
input_correlation_mat_bg = paste0(hg38_gene, '.', mm10_gene_bg_i, '.cor.heatmap.png.cor.mat.txt')
if (file.exists(input_correlation_mat_bg)) {
# read in correlation matrix
cor_mat_bg = read.table(input_correlation_mat_bg, header = F)
if (dim(cor_mat_bg)[1]!=dim(cor_mat)[1]){
    next
}
print(dim(cor_mat_bg))

# smooth cor_mat matrix along the columns for each row
cor_mat_bg_smooth = smooth_cor_mat(cor_mat_bg, spar = smooth_spar)

# negative control Cd4 gene
cor_mat_bg_zcol = matrix(0, nrow=dim(cor_mat)[1], ncol=dim(cor_mat)[2])
non0_col = c()
for (j in 1:dim(cor_mat)[2]) {
  # get cor_mat high signals based zpfdr
  if (cor_mat_bg_r_smooth_mean[j]!=0){
    cor_mat_bg_zcol[,j] = (cor_mat_bg_smooth[,j] - cor_mat_bg_r_smooth_mean[j]) / cor_mat_bg_r_smooth_sd[j]
  }else {
    cor_mat_bg_zcol[,j] = (cor_mat_bg_smooth[,j] - mean(cor_mat_bg_smooth[,j])) / sd(cor_mat_bg_smooth[,j])
  }
}

cor_mat_bg_zpfdr = apply(cor_mat_bg_zcol, 2, function(x) p.adjust(pnorm(x, lower.tail=F), method = 'fdr'))
FP_n = c(FP_n, sum(cor_mat_bg_zpfdr<fdr_thresh))
all_bgn = c(all_bgn, dim(cor_mat_bg)[1]*dim(cor_mat_bg)[2])
print(sum(cor_mat_bg_zpfdr<fdr_thresh))

# naive FDR
cor_mat_bg_local_raw = get_2r_fdr(cor_mat_bg_smooth, fdr_thresh)
FP_n_naive = c(FP_n_naive, sum(cor_mat_bg_local_raw!=0))
}}

# get FDR random picked bg genes adjusted
TP = sum(cor_mat_zpfdr<fdr_thresh)
FP_ave = mean(FP_n)
FRD = FP_ave/(TP+FP_ave)
print(FRD)

# get FDR naive
TP = sum(cor_mat_zpfdr<fdr_thresh)
FP_naive_ave = mean(FP_n_naive)
FRD_naive = FP_naive_ave/(TP+FP_naive_ave)
print(FRD_naive)


