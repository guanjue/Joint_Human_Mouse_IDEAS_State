setwd('/Users/guanjuexiang/Documents/projects/analysis/0813_human_mouse_state_compare_heatmap')

library(LSD)

hg38_gene = 'GATA1'
mm10_gene = 'Gata1'
mm10_gene_bg = 'Cd4'
input_correlation_mat = paste0('GATA1.Gata1.cor.heatmap.png.cor.mat.RAW_data_for_plot.txt')
input_correlation_mat_bg = paste0(hg38_gene, '.', mm10_gene_bg, '.cor.heatmap.png.cor.mat.txt')
input_sequence_similarity_mat = paste0('GATA1_hg38_Gata1_mm10_DNA_sequence_similarity.txt')

fdr_thresh = 1e-1

# read in correlation matrix
cor_mat = read.table(input_correlation_mat, header = F)
# read in correlation matrix
cor_mat_bg = read.table(input_correlation_mat_bg, header = F)

# plot heatscatter of cor_mat
cor_mat_heatscatter_x = c()
cor_mat_heatscatter_y = c()
for (i in 1:nrow(cor_mat)) {
    cor_mat_heatscatter_x = c(cor_mat_heatscatter_x, 1:ncol(cor_mat))
    cor_mat_heatscatter_y = c(cor_mat_heatscatter_y, as.numeric(cor_mat[i,]))
}

# plot cor_mat along the columns
png(paste0(hg38_gene, '.', mm10_gene, '.cor_column.cor.png'), width=1000, height=400)
heatscatter(cor_mat_heatscatter_x, cor_mat_heatscatter_y, pch = 20, cex = 1, xlab = 'Column', ylab = 'Correlation', ylim=c(min(0),max((1))))
#abline(h = min(cor_mat[cor_mat_zpfdr<0.1]), col = 'red')
dev.off()

smooth_cor_mat = function(cor_mat, spar = 0.1) {
  cor_mat_smooth_col = t(apply(cor_mat, 1, function(x) smooth.spline(1:ncol(cor_mat), x, spar = spar)$y))
  cor_mat_smooth_row = apply(cor_mat, 2, function(x) smooth.spline(1:nrow(cor_mat), x, spar = spar)$y)
  cor_mat_smooth = (cor_mat_smooth_col + cor_mat_smooth_row)/2
  return(cor_mat_smooth)
}

cor_mat_smooth = smooth_cor_mat(cor_mat, spar = 0.1)
cor_mat_bg_smooth = smooth_cor_mat(cor_mat_bg, spar = 0.1)

# read random_pick_bg_genes
random_picked_bg_genes = 'randomly_pick_bg_genes/mm10.gene.all.txt'
random_picked_bg_genes = read.table(random_picked_bg_genes, header = F)

# initial a list to store cor_mat_bg_r_smooth
cor_mat_bg_r_smooth_mean = c()
cor_mat_bg_r_smooth_sd = c()
# read cor_mat for random picked bg genes
cor_mat_bg_r_smooth_mat = matrix(-1000, nrow=dim(random_picked_bg_genes)[1]*dim(cor_mat_bg_r)[1], ncol=dim(cor_mat_bg_r)[2])
cor_mat_bg_r_smooth_mat_start = rep(1, dim(cor_mat_bg)[2])

for (i in 1:nrow(random_picked_bg_genes)) {
  mm10_gene_bg = random_picked_bg_genes[i,1]
  print(mm10_gene_bg)
  #mm10_gene_bg = 'Cd4'
  input_correlation_mat_bg = paste0(hg38_gene, '.', mm10_gene_bg, '.cor.heatmap.png.cor.mat.txt')
  if (file.exists(input_correlation_mat_bg)) {
  cor_mat_bg_r = read.table(input_correlation_mat_bg, header = F)
  cor_mat_bg_r_smooth = smooth_cor_mat(cor_mat_bg, spar = 0.1)
  # get fdr top genes
  for (j in 1:dim(cor_mat_bg_r_smooth)[2]){
    cor_mat_bg_r_smooth_fdr_j = p.adjust(pnorm((cor_mat_bg_r_smooth[,j]-mean(cor_mat_bg_r_smooth[,j]))/sd(cor_mat_bg_r_smooth[,j]), lower.tail=F), method = 'fdr')
    # get cor_mat high signals based zpfdr
    #print(sum(cor_mat_bg_r_smooth_fdr_j<fdr_thresh))
    if (sum(cor_mat_bg_r_smooth_fdr_j<fdr_thresh)>1){
      cor_mat_bg_r_smooth_mat[cor_mat_bg_r_smooth_mat_start[j]:(cor_mat_bg_r_smooth_mat_start[j]+sum(cor_mat_bg_r_smooth_fdr_j<fdr_thresh)-1),j] = cor_mat_bg_r_smooth[cor_mat_bg_r_smooth_fdr_j<fdr_thresh,j]
      cor_mat_bg_r_smooth_mat_start[j] = cor_mat_bg_r_smooth_mat_start[j] + sum(cor_mat_bg_r_smooth_fdr_j<fdr_thresh) 
    } 
  }
  print(sum(cor_mat_bg_r_smooth_mat_start!=1))
}
}

cor_mat_bg_r_smooth_mean = apply(cor_mat_bg_r_smooth_mat, 2, function(x) mean(x[x!=-1000]))
cor_mat_bg_r_smooth_sd = apply(cor_mat_bg_r_smooth_mat, 2, function(x) sd(x[x!=-1000]))
cor_mat_bg_r_smooth_mean[is.na(cor_mat_bg_r_smooth_mean)] = 0
cor_mat_bg_r_smooth_sd[is.na(cor_mat_bg_r_smooth_sd)] = 0

# FDR along the columns
fdr_thresh = 1e-1
cor_mat_zcol = matrix(0, nrow=dim(cor_mat)[1], ncol=dim(cor_mat)[2])
non0_col = c()
for (j in 1:dim(cor_mat)[2]) {
  # get cor_mat high signals based zpfdr
  print(cor_mat_bg_r_smooth_mean[j]!=0)
  if (cor_mat_bg_r_smooth_mean[j]!=0){
    cor_mat_zcol[,j] = (cor_mat_smooth[,j] - cor_mat_bg_r_smooth_mean[j]) / cor_mat_bg_r_smooth_sd[j]
  }else {
    cor_mat_zcol[,j] = (cor_mat_smooth[,j] - mean(cor_mat_smooth[,j])) / sd(cor_mat_smooth[,j])
  }
}

cor_mat_zpfdr = apply(cor_mat_zcol, 2, function(x) p.adjust(pnorm(x, lower.tail=F), method = 'fdr'))
#cor_mat_zpfdr = apply(cor_mat_zcol, 2, function(x) (pnorm(x, lower.tail=F)))
cor_mat_local = cor_mat_smooth * (cor_mat_zpfdr<fdr_thresh)
# write cor_mat_local
write.table(cor_mat_local, paste0(hg38_gene, '.', mm10_gene, '.cor.heatmap.png.cor.mat.local.txt'), row.names = F, col.names = F, quote = F, sep = '\t')

# bg
fdr_thresh = 1e-1
cor_mat_bg_zcol = matrix(0, nrow=dim(cor_mat)[1], ncol=dim(cor_mat)[2])
non0_col = c()
for (j in 1:dim(cor_mat)[2]) {
  # get cor_mat high signals based zpfdr
  print(cor_mat_bg_r_smooth_mean[j]!=0)
  if (cor_mat_bg_r_smooth_mean[j]!=0){
    cor_mat_bg_zcol[,j] = (cor_mat_bg_smooth[,j] - cor_mat_bg_r_smooth_mean[j]) / cor_mat_bg_r_smooth_sd[j]
  }else {
    cor_mat_bg_zcol[,j] = (cor_mat_bg_smooth[,j] - mean(cor_mat_bg_smooth[,j])) / sd(cor_mat_bg_smooth[,j])
  }
}

cor_mat_bg_zpfdr = apply(cor_mat_bg_zcol, 2, function(x) p.adjust(pnorm(x, lower.tail=F), method = 'fdr'))
#cor_mat_zpfdr = apply(cor_mat_zcol, 2, function(x) (pnorm(x, lower.tail=F)))
cor_mat_bg_local = cor_mat_bg_smooth * (cor_mat_bg_zpfdr<fdr_thresh)
# write cor_mat_local
write.table(cor_mat_bg_local, paste0(hg38_gene, '.', mm10_gene, '.cor.heatmap.png.cor.mat.local.bg.txt'), row.names = F, col.names = F, quote = F, sep = '\t')








# FDR along the columns
fdr_thresh = 1e-1
cor_mat_zcol = matrix(0, nrow=dim(cor_mat)[1], ncol=dim(cor_mat)[2])
non0_col = c()
for (j in 1:dim(cor_mat)[2]) {
  # get cor_mat_bg zpfdr
  cor_mat_bg_zpfdr = p.adjust(pnorm((cor_mat_bg_smooth[,j] - mean(cor_mat_bg_smooth[,j])) / sd(cor_mat_bg_smooth[,j]), lower.tail=F), method = 'fdr')
  cor_mat_bg_zp = (pnorm((cor_mat_bg[,j] - mean(cor_mat_bg[,j])) / sd(cor_mat_bg[,j]), lower.tail=F))
  # get cor_mat high signals based zpfdr
  print(sum(cor_mat_bg_zpfdr<fdr_thresh))
  if (sum(cor_mat_bg_zpfdr<fdr_thresh)>1){
    cor_mat_zcol[,j] = (cor_mat_smooth[,j] - mean(cor_mat_bg_smooth[cor_mat_bg_zpfdr<fdr_thresh,j])) / sd(cor_mat_bg_smooth[cor_mat_bg_zpfdr<fdr_thresh,j])
    non0_col = c(non0_col, sum(cor_mat_bg_zpfdr<fdr_thresh))
  }else {
    cor_mat_zcol[,j] = (cor_mat_smooth[,j] - mean(cor_mat_bg_smooth[,j])) / sd(cor_mat_bg_smooth[,j])
  }
}

cor_mat_zpfdr = apply(cor_mat_zcol, 2, function(x) p.adjust(pnorm(x, lower.tail=F), method = 'fdr'))
#cor_mat_zpfdr = apply(cor_mat_zcol, 2, function(x) (pnorm(x, lower.tail=F)))
cor_mat_local = cor_mat * (cor_mat_zpfdr<fdr_thresh)
# write cor_mat_local
write.table(cor_mat_local, paste0(hg38_gene, '.', mm10_gene, '.cor.heatmap.png.cor.mat.local.txt'), row.names = F, col.names = F, quote = F, sep = '\t')



png(paste0(hg38_gene, '.', mm10_gene, '.cor_column.cor.r.png'), width=1000, height=400)
plot(1:ncol(cor_mat), (cor_mat_smooth[1,]), pch = 20, cex = 1, xlab = 'Column', ylab = 'Correlation', ylim=c(min(0),max((1))))
for (i in 2:nrow(cor_mat)) {
  points(1:ncol(cor_mat), (cor_mat_smooth[i,]), pch = 20, cex = 1)
}
for (i in 1:nrow(cor_mat)) {
  points((1:ncol(cor_mat))[cor_mat_local[i,]!=0], (cor_mat_local[i,][cor_mat_local[i,]!=0]), pch = 20, cex = 1, col='red')
}
dev.off()

png(paste0(hg38_gene, '.', mm10_gene, '.cor_column.cor.r.bg.png'), width=1000, height=400)
plot(1:ncol(cor_mat_bg), (cor_mat_bg_smooth[1,]), pch = 20, cex = 1, xlab = 'Column', ylab = 'Correlation', ylim=c(min(0),max((1))))
for (i in 2:nrow(cor_mat_bg)) {
  points(1:ncol(cor_mat_bg), (cor_mat_bg_smooth[i,]), pch = 20, cex = 1)
}
for (i in 1:nrow(cor_mat_bg)) {
  points((1:ncol(cor_mat_bg))[cor_mat_bg_local[i,]!=0], (cor_mat_bg_local[i,][cor_mat_bg_local[i,]!=0]), pch = 20, cex = 1, col='red')
}
dev.off()












# plot cor_mat_zpfdr along the columns
png(paste0(hg38_gene, '.', mm10_gene, '.cor_column.zpfdr.png'), width=1000, height=500)
plot(1:ncol(cor_mat), -log10(cor_mat_zpfdr[1,]), pch = 20, cex = 1, xlab = 'Column', ylab = 'FDR adjusted Z score based p-value', ylim=c(0,max(-log10(cor_mat_zpfdr))))
for (i in 2:nrow(cor_mat)) {
  points(1:ncol(cor_mat), -log10(cor_mat_zpfdr[i,]), pch = 20, cex = 1)
}
abline(h = -log10(0.1), col = 'red')
dev.off()

# smooth cor_mat matrix along the columns for each row
cor_mat_smooth = t(apply(cor_mat, 1, function(x) smooth.spline(1:ncol(cor_mat), x, spar = 0.1)$y))
cor_mat_smooth_zcol = apply(cor_mat_smooth, 2, function(x) (x - mean(x)) / sd(x))
cor_mat_smooth_zpfdr = apply(cor_mat_smooth_zcol, 2, function(x) p.adjust(pnorm(x, lower.tail=F), method = 'fdr'))
cor_zpfdr_thresh = max(cor_mat_smooth_zpfdr[cor_mat_smooth >= cor_thresh])
fdr_thresh = 1e-1
cor_mat_smooth_local = cor_mat_smooth * (cor_mat_zpfdr<fdr_thresh)

# smooth cor_mat matrix along the columns for each row
png(paste0(hg38_gene, '.', mm10_gene, '.cor_column.cor_sm.png'), width=1000, height=500)
plot(1:ncol(cor_mat), (cor_mat_smooth[1,]), pch = 20, cex = 1, xlab = 'Column', ylab = 'Correlation', ylim=c(min(0),max((cor_mat_smooth))))
for (i in 2:nrow(cor_mat)) {
  points(1:ncol(cor_mat), (cor_mat_smooth[i,]), pch = 20, cex = 1)
}
for (i in 1:nrow(cor_mat)) {
  points(1:ncol(cor_mat), (cor_mat_smooth_local[i,]), pch = 20, cex = 1, col='red')
}
#abline(h = min(cor_mat_smooth[cor_mat_smooth_zpfdr<1e-3]), col = 'red')
dev.off()

png(paste0(hg38_gene, '.', mm10_gene, '.cor_column.zpfdr.sm.png'), width=1000, height=500)
plot(1:ncol(cor_mat), -log10(cor_mat_smooth_zpfdr[1,]), pch = 20, cex = 1, xlab = 'Column', ylab = 'FDR adjusted Z score based p-value', ylim=c(0,max(-log10(cor_mat_smooth_zpfdr))))
for (i in 2:nrow(cor_mat)) {
  points(1:ncol(cor_mat), -log10(cor_mat_smooth_zpfdr[i,]), pch = 20, cex = 1)
}
abline(h = -log10(0.1), col = 'red')
dev.off()



# plot scatter plot of correlation vs sequence similarity
png(paste0(hg38_gene, '.', mm10_gene, '.cor_vs_seq.png'))
heatscatter(as.numeric(as.matrix(seq_mat)), as.numeric(as.matrix(cor_mat)), xlab = 'Sequence similarity', ylab = 'Correlation', pch = 20, cex = 0.5)
dev.off()

png(paste0(hg38_gene, '.', mm10_gene, '.cor_vs_seq.max.png'))
heatscatter(as.numeric(c(apply(seq_mat,2,max), apply(seq_mat,1,max))), as.numeric(c(apply(cor_mat,2,max), apply(cor_mat,1,max))), xlab = 'Sequence similarity', ylab = 'Correlation', pch = 20, cex = 0.5)
dev.off()

# pheatmap white black
png(paste0(hg38_gene, '.', mm10_gene, '.cor_vs_seq.pheatmap.png'))
pheatmap(cor_mat, color = colorRampPalette(c('white', 'black'))(100), cluster_rows = F, cluster_cols = F, fontsize = 8, border_color = NA, cellwidth = 1, cellheight = 1)
dev.off()

