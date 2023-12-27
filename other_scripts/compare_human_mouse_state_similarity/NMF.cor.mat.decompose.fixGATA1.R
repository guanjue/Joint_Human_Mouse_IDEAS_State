setwd('/Users/guanjuexiang/Documents/projects/analysis/0813_human_mouse_state_compare_heatmap')

# conda activate shiny
library(pheatmap)
#library(ComplexHeatmap)
#library(InteractiveComplexHeatmap)
#library(shiny)
#library(circlize)

# define NMF function
nmf <- function(V, k, max_iter = 1000, tol = 1e-4) {
  # Initialize W and H with random values
  n <- dim(V)[1]
  m <- dim(V)[2]
  W <- matrix(runif(n * k), nrow = n)
  H <- matrix(runif(k * m), nrow = k)
  for (i in 1:max_iter) {
    # Update H
    WH <- W %*% H
    H <- H * ((t(W) %*% (V / (WH + 1e-10))) / (t(W) %*% matrix(1, nrow = n, ncol = m)))
    # Update W
    WH <- W %*% H
    W <- W * (((V / (WH + 1e-10)) %*% t(H)) / (matrix(1, nrow = n, ncol = m) %*% t(H)))
    # Normalize W and H
    if (i < round(max_iter*0.5)) {
      W <- W / max(W)
      H <- H / max(H)
    }
    # Check convergence
    WH <- W %*% H
    obj <- sum(V * log(1 / (WH + 1e-10)) + WH)
    if (i > 1 && abs(obj - prev_obj) < tol) {
      break
    }
    prev_obj <- obj
  }  
  output = list(W = W, H = H)
  return(output)
}

nmf_fix_H <- function(V, H, max_iter = 1000, tol = 1e-4) {
  # Initialize W and H with random values
  n <- dim(V)[1]
  m <- dim(V)[2]
  k <- dim(H)[1]
  W <- matrix(runif(n * k), nrow = n)
  H <- H #matrix(runif(k * m), nrow = k)
  for (i in 1:max_iter) {
    # Update H
    WH <- W %*% H
    #H <- H * ((t(W) %*% (V / (WH + 1e-10))) / (t(W) %*% matrix(1, nrow = n, ncol = m)))
    # Update W
    WH <- W %*% H
    W <- W * (((V / (WH + 1e-10)) %*% t(H)) / (matrix(1, nrow = n, ncol = m) %*% t(H)))
    # Normalize W and H
    if (i < round(max_iter*0.5)) {
      W <- W / max(W)
      #H <- H / max(H)
    }
    # Check convergence
    WH <- W %*% H
    obj <- sum(V * log(1 / (WH + 1e-10)) + WH)
    if (i > 1 && abs(obj - prev_obj) < tol) {
      break
    }
    prev_obj <- obj
  }  
  output = list(W = W, H = H)
  return(output)
}


zp2 = function(x, thresh) {
    # 1st round z score
    z = (x - mean(x)) / sd(x)
    # get z score one-side p value
    zp = 1 - pnorm(z)
    # 2nd round z score
    z2 = (x - mean(x[zp>thresh])) / sd(x[zp>thresh])
    # get z score one-side p value
    zp2 = 1 - pnorm(z2)
    return(zp2)
}

zp = function(x, thresh) {
    # 1st round z score
    z = (x - mean(x)) / sd(x)
    # get z score one-side p value
    zp = 1 - pnorm(z)
    return(zp)
}

zp2_bg = function(x, xbg) {
    # 1st round z score
    z = (x - mean(xbg)) / sd(xbg)
    # get z score one-side p value
    zp = 1 - pnorm(z)
    return(zp)
}


zp2_sm = function(x, thresh, span) {
    x_sm = loess(x ~ c(1:length(x)), span=span)$fitted
    x_sm[x_sm<=0] = min(x_sm[x_sm>0])
    x_sm_zp = zp2(x_sm, thresh)
    x_sm_zp_FDR = p.adjust(x_sm_zp, method = 'fdr')
    return(x_sm_zp_FDR)
}


zp_sm = function(x, thresh, span) {
    cor_factor = x
    cor_factor_sm = loess(cor_factor ~ c(1:length(x)), span=span)$fitted
    cor_factor_sm[cor_factor_sm<=0] = min(cor_factor_sm[cor_factor_sm>0])
    cor_factor_sm_zp = zp(cor_factor_sm, thresh)
    cor_factor_sm_zp_FDR = p.adjust(cor_factor_sm_zp, method = 'fdr')
    return(cor_factor_sm_zp_FDR)
}

zp2_sm_bg = function(x, xbg, span) {
    # read x
    x_sm = loess(x ~ c(1:length(x)), span=span)$fitted
    x_sm[x_sm<=0] = min(x_sm[x_sm>0])
    x_sm_zp = zp2_bg(x_sm, xbg)
    return(x_sm_zp)
}

zp2_sm_bg_fdr = function(x, xbg, span) {
    # read x
    x_sm = loess(x ~ c(1:length(x)), span=span)$fitted
    x_sm[x_sm<=0] = min(x_sm[x_sm>0])
    x_sm_zp = zp2_bg(x_sm, xbg)
    # get FDR
    x_sm_zp_FDR = p.adjust(x_sm_zp, method = 'fdr')
    return(x_sm_zp_FDR)
}

get_FDR_p_thresh = function(nmf_result_W_j, nmf_result_W_mat_long_j, start, end, bin, FDR_thresh){
    nmf_result_W_sm_j = loess(nmf_result_W_j ~ c(1:length(nmf_result_W_j)), span=0.01)$fitted
    nmf_result_W_mat_long_sm_j = loess(nmf_result_W_mat_long_j ~ c(1:length(nmf_result_W_mat_long_j)), span=0.01)$fitted
    for (p_i in seq(start, end, by=bin)){
        if (p_i%%0.01==0){print(p_i)}
        #print((mean(nmf_result_W_sm_j<p_i))/(mean(nmf_result_W_mat_long_sm_j<p_i)))
        FDR_p_i = mean(nmf_result_W_mat_long_sm_j<p_i)/(mean(nmf_result_W_mat_long_sm_j<p_i)+mean(nmf_result_W_sm_j<p_i))
        if (FDR_p_i<FDR_thresh) {
            print(c(p_i, FDR_p_i))
            FDR_p_thresh = p_i
            break
        }
    }
    return(FDR_p_thresh)
}


# output plot folder
hg38_gene = 'GATA1'
mm10_gene = 'Gata1'
mm10_gene_bg = 'Gata1'
mm10_gene_bg = 'Cd4'
mm10_gene_bg = 'Rps19'
mm10_gene_bg = 'Slc4a1'

mm10_gene_bg = 'Slc4a1'
output_folder = paste0('NMF_reconstruction_fix_GATA_', hg38_gene, '_', mm10_gene_bg)
input_correlation_mat = paste0(hg38_gene, '.', mm10_gene, '.cor.heatmap.png.cor.mat.txt')
input_correlation_mat_bg = paste0(hg38_gene, '.', mm10_gene_bg, '.cor.heatmap.png.cor.mat.txt')
hg38_coord = paste0('hg38.gene.', hg38_gene, '.matched_ct.state.bed')
mm10_coord = paste0('mm10.gene.', mm10_gene, '.matched_ct.state.bed')
random_picked_bg_genes_file = 'randomly_pick_bg_genes/mm10.gene.all.txt'

fdr_thresh = 1e-1


# read correlation matrix
d_mk_cor = as.matrix(read.table(input_correlation_mat, header=F, sep='\t', comment.char='~'))
d_mk_cor_mat = as.matrix(d_mk_cor)
d_mk_cor_bg = as.matrix(read.table(input_correlation_mat_bg, header=F, sep='\t', comment.char='~'))
d_mk_cor_mat_bg = as.matrix(d_mk_cor_bg)

# read coordinates
bed_hg38 = read.table(hg38_coord, header=F, sep='\t', comment.char='~')
bed_mm10 = read.table(mm10_coord, header=F, sep='\t', comment.char='~')
# add rownames and colnames
colnames(d_mk_cor_mat) = apply(bed_hg38, 1, function(x) paste0('H_', x[1], '_', x[2], '_', x[3]))
rownames(d_mk_cor_mat) = apply(bed_mm10, 1, function(x) paste0('M_', x[1], '_', x[2], '_', x[3]))
colnames(d_mk_cor_mat_bg) = apply(bed_hg38, 1, function(x) paste0('H_', x[1], '_', x[2], '_', x[3]))
rownames(d_mk_cor_mat_bg) = apply(bed_mm10, 1, function(x) paste0('M_', x[1], '_', x[2], '_', x[3]))

# get positive correlation matrix
d_mk_cor_mat_pos = d_mk_cor_mat
d_mk_cor_mat_pos[d_mk_cor_mat_pos<0] = 0
d_mk_cor_mat_bg_pos = d_mk_cor_mat_bg
d_mk_cor_mat_bg_pos[d_mk_cor_mat_bg_pos<0] = 0

# save positive correlation matrix

# mkdir
dir.create(output_folder)

#  NMF decomposition
set.seed(2019)
n_comp = 6
nmf_result_0 <- nmf(d_mk_cor_mat_pos, k = n_comp)

# Extract the independent components (ICs) and the mixing matrix
nmf_result_0_W <- nmf_result_0$W
nmf_result_0_H <- nmf_result_0$H

##############################################################################################################

# plot cross factor correlation matrix
# get top quantile of NMFs
nmf_result_0_W_binary_top_quantile = apply(nmf_result_0_W, 2, function(x) x > quantile(x, 0.9))*1
nmf_result_0_H_binary_top_quantile = apply(t(nmf_result_0_H), 2, function(x) x > quantile(x, 0.9))*1

# get mean correlation between different NMF factors
cross_factor_cor_mat = matrix(0, nrow = dim(nmf_result_0_W)[2], ncol = dim(nmf_result_0_W)[2])
#
for (hg38_factor_j in 1:dim(nmf_result_0_H)[1]){
    for (mm10_factor_i in 1:dim(nmf_result_0_W)[2]){
        # get FDR_hg38_j and FDR_mm10_i
        FDR_hg38_j = nmf_result_0_H_binary_top_quantile[,hg38_factor_j]
        FDR_mm10_i = nmf_result_0_W_binary_top_quantile[,mm10_factor_i]
        # get positive correlation matrix
        d_mk_cor_mat_pos_cols = d_mk_cor_mat_pos[,FDR_hg38_j!=0]
        d_mk_cor_mat_pos_cols_rows = d_mk_cor_mat_pos_cols[FDR_mm10_i!=0,]
        # get mean correlation
        cross_factor_cor_mat[mm10_factor_i,hg38_factor_j] = mean(d_mk_cor_mat_pos_cols_rows)
    }
}

# plot cross_factor_cor_mat
cross_factor_cor_mat[is.na(cross_factor_cor_mat)] = 0
colnames(cross_factor_cor_mat) = paste0(hg38_gene, '_', 1:dim(nmf_result_0_W)[2])
rownames(cross_factor_cor_mat) = paste0(mm10_gene, '_', 1:dim(nmf_result_0_W)[2])
# my colorbar white to red
my_colorbar = colorRampPalette(c('white','red'))(n = 100)
pdf(paste0(hg38_gene, '.', mm10_gene, '.cross_factor_cor_mat.pdf'), width = 10, height = 10)
pheatmap(cross_factor_cor_mat, cluster_rows=F, cluster_cols=F, color = my_colorbar, legend=T, show_rownames=T, show_colnames=T, cex=1.5)
dev.off()

# write cross_factor_cor_mat
write.table(cross_factor_cor_mat, paste0(hg38_gene, '.', mm10_gene, '.cross_factor_cor_mat.txt'), sep='\t', quote=F, row.names=T, col.names=T)

##############################################################################################################

# nmf fix H
set.seed(2019)
nmf_result <- nmf_fix_H(d_mk_cor_mat_bg_pos, H = nmf_result_0_H)
# Extract the independent components (ICs) and the mixing matrix
nmf_result_W <- nmf_result$W
nmf_result_H <- nmf_result$H


# cluster factor matrix
nmf_W_H = rbind(nmf_result_0_W, t(nmf_result_0_H))
nmf_W_H_order_tree = hclust(dist(t(nmf_W_H)), method = 'average')
nmf_W_H_order = nmf_W_H_order_tree$order

# plot NMFs_mm10 heatmap
#my_colorbar=colorRampPalette(c('black', 'white'))(n = 100)
png(paste0(output_folder, '/', hg38_gene, '.', mm10_gene, '.cor.heatmap.png.cor.mat.txt.NMFs_mm10.heatmap.png'), width = 300, height = 1000)
pheatmap(nmf_result_W[,], cluster_rows=F, cluster_cols=F, width = 2, height = 6, legend=F, show_rownames=F, show_colnames=F)
dev.off()
# plot NMFs_hg38 heatmap
png(paste0(output_folder, '/', hg38_gene, '.', mm10_gene, '.cor.heatmap.png.cor.mat.txt.NMFs_hg38.heatmap.png'), width = 300, height = 1000)
pheatmap(t(nmf_result_H)[,], cluster_rows=F, cluster_cols=F, width = 2, height = 6, legend=F, show_rownames=F, show_colnames=F)
dev.off()

# naive zp plot for each NMF factor
nmf_result_0_W_binary = nmf_result_0_W
for (j in 1:dim(nmf_result_0_W)[2]){
    nmf_result_0_W_sm_zpfdr = (zp2_sm(nmf_result_0_W[,j], 0.01, 0.03) < fdr_thresh) *1
    nmf_result_0_W_binary[,j] = nmf_result_0_W_sm_zpfdr
}
nmf_result_W_binary = nmf_result_W
for (j in 1:dim(nmf_result_W)[2]){
    nmf_result_W_sm_zpfdr = (zp2_sm(nmf_result_W[,j], 0.01, 0.03) < fdr_thresh) *1
    nmf_result_W_binary[,j] = nmf_result_W_sm_zpfdr
}

# plot NMFs_mm10 binary heatmap
my_colorbar_binary=colorRampPalette(c('white','black'))(n = 100)
png(paste0(output_folder, '/', hg38_gene, '.', mm10_gene, '.cor.heatmap.png.cor.mat.txt.NMFs_mm10.binary_heatmap.png'), width = 300, height = 1000)
pheatmap(nmf_result_0_W_binary, cluster_rows=F, cluster_cols=F, width = 2, height = 6, legend=F, show_rownames=F, show_colnames=F, color=my_colorbar_binary)
dev.off()
#
png(paste0(output_folder, '/', hg38_gene, '.', mm10_gene_bg, '.cor.heatmap.png.cor.mat.txt.NMFs_mm10.binary_heatmap.png'), width = 300, height = 1000)
pheatmap(nmf_result_W_binary, cluster_rows=F, cluster_cols=F, width = 2, height = 6, legend=F, show_rownames=F, show_colnames=F, color=my_colorbar_binary)
dev.off()


# FDR of each NMF factor for foreground
nmf_result_0_W_binary_FDR = (colSums(nmf_result_W_binary)) / (colSums(nmf_result_W_binary)+colSums(nmf_result_0_W_binary))
nmf_result_0_W_binary_FDR

##############################################################################################################
# read random_pick_bg_genes
random_picked_bg_genes = read.table(random_picked_bg_genes_file, header = F)

# read cor_mat for random picked bg genes
nmf_result_W_mat_long = c()

for (i in 1:nrow(random_picked_bg_genes)) {
    # read gene name
    mm10_gene_bg_i = random_picked_bg_genes[i,1]
    print(mm10_gene_bg_i)
    # get mm10_gene_bg file
    input_correlation_mat_bg_i = paste0(hg38_gene, '.', mm10_gene_bg_i, '.cor.heatmap.png.cor.mat.txt')

    if (file.exists(input_correlation_mat_bg_i)) {
    # read correlation matrix
    d_mk_cor_bg_i = as.matrix(read.table(input_correlation_mat_bg_i, header=F, sep='\t', comment.char='~'))
    d_mk_cor_mat_bg_i = as.matrix(d_mk_cor_bg_i)
    #
    if ((dim(d_mk_cor_mat_bg_i)[1] != dim(bed_mm10)[1]) | (dim(d_mk_cor_mat_bg_i)[2] != dim(bed_hg38)[1])) {
        next
    }
    # add rownames and colnames
    colnames(d_mk_cor_mat_bg_i) = apply(bed_hg38, 1, function(x) paste0('H_', x[1], '_', x[2], '_', x[3]))
    rownames(d_mk_cor_mat_bg_i) = apply(bed_mm10, 1, function(x) paste0('M_', x[1], '_', x[2], '_', x[3]))

    # get positive correlation matrix
    d_mk_cor_mat_bg_pos_i = d_mk_cor_mat_bg_i
    d_mk_cor_mat_bg_pos_i[d_mk_cor_mat_bg_pos_i<0] = 0

    # nmf fix H
    set.seed(2019)
    nmf_result_i <- nmf_fix_H(d_mk_cor_mat_bg_pos_i, H = nmf_result_0_H)
    # Extract the independent components (ICs) and the mixing matrix
    nmf_result_W_i <- nmf_result_i$W

    # add to nmf_result_W_mat_long
    nmf_result_W_mat_long = rbind(nmf_result_W_mat_long, nmf_result_W_i)
    # get significant high index
    #nmf_result_W_i_binary = matrix(0, nrow = dim(nmf_result_W_i)[1], ncol = dim(nmf_result_W_i)[2])
    #for (j in 1:dim(nmf_result_0_W)[2]){
    #    nmf_result_W_i_sm_zpfdr = (zp2_sm(nmf_result_W_i[,j], 0.01, 0.03) < 0.1) *1
    #    if (sum(nmf_result_W_i_sm_zpfdr) != 0) {
    #        nmf_result_W_i_binary[,j] = nmf_result_W_i_sm_zpfdr
    #    }
    #}

    # add to nmf_result_W_mat_long
    #nmf_result_W_mat_long = rbind(nmf_result_W_mat_long, nmf_result_W_i*nmf_result_W_i_binary)
    #print(apply(nmf_result_W_i*nmf_result_W_i_binary,2,function(x) mean(x[x!=0])))
    #print(apply(nmf_result_W_mat_long,2,function(x) mean(x[x!=0])))
}}
##############################################################################################################

# get binary matrix based on random selected bg genes for foreground gene
nmf_result_0_W_binary = nmf_result_0_W
FDR_p_thresh_0_vec = c()
for (j in 1:dim(nmf_result_0_W)[2]){
    #nmf_result_0_W_sm_zpfdr = (zp2_sm_bg_fdr(nmf_result_0_W[,j], nmf_result_W_mat_long[,j], 0.03) < fdr_thresh) *1
    nmf_result_0_W_j = nmf_result_0_W[,j]
    nmf_result_W_mat_long_j = nmf_result_W_mat_long[,j]
    # get FDR_p_thresh
    FDR_p_thresh_0_j = get_FDR_p_thresh(nmf_result_0_W_j, nmf_result_W_mat_long_j, 0.05, 1e-10, -0.0001, fdr_thresh)
    FDR_p_thresh_0_vec = c(FDR_p_thresh_0_vec, FDR_p_thresh_0_j)
    # get binary matrix
    nmf_result_0_W_sm_zpfdr = (zp2_sm_bg(nmf_result_0_W_j, nmf_result_W_mat_long_j, 0.03) < FDR_p_thresh_0_j) *1
    nmf_result_0_W_binary[,j] = nmf_result_0_W_sm_zpfdr
}

# get binary matrix based on random selected bg genes for control gene
nmf_result_W_binary = nmf_result_W
for (j in 1:dim(nmf_result_W)[2]){
    #nmf_result_W_sm_zpfdr = (zp2_sm_bg_fdr(nmf_result_W[,j], nmf_result_W_mat_long[,j], 0.03) < fdr_thresh) *1
    nmf_result_W_sm_zpfdr = (zp2_sm_bg(nmf_result_W[,j], nmf_result_W_mat_long[,j], 0.03) < FDR_p_thresh_0_vec[j]) *1
    nmf_result_W_binary[,j] = nmf_result_W_sm_zpfdr
}

# plot NMFs_mm10 binary heatmap
my_colorbar_binary=colorRampPalette(c('white','black'))(n = 100)
png(paste0(output_folder, '/', hg38_gene, '.', mm10_gene, '.cor.heatmap.png.cor.mat.txt.NMFs_mm10.binary_bgadj_heatmap.png'), width = 300, height = 1000)
pheatmap(nmf_result_0_W_binary, cluster_rows=F, cluster_cols=F, width = 2, height = 6, legend=F, show_rownames=F, show_colnames=F, color=my_colorbar_binary)
dev.off()

my_colorbar_binary=colorRampPalette(c('white','black'))(n = 100)
png(paste0(output_folder, '/', hg38_gene, '.', mm10_gene_bg, '.cor.heatmap.png.cor.mat.txt.NMFs_mm10.binary_bgadj_heatmap.png'), width = 300, height = 1000)
pheatmap(nmf_result_W_binary, cluster_rows=F, cluster_cols=F, width = 2, height = 6, legend=F, show_rownames=F, show_colnames=F, color=my_colorbar_binary)
dev.off()

# FDR of each NMF factor for foreground
nmf_result_0_W_binary_FDR = (colSums(nmf_result_W_binary)) / (colSums(nmf_result_W_binary)+colSums(nmf_result_0_W_binary))
nmf_result_0_W_binary_FDR
mm10_gene_bg

##############################################################################################################

Cd4_naive_FDR = c(NA, 0.5964912, 0.3625000, 0.5882353, NA, 0.3364486)
Rps19_naive_FDR = c(NA, 0.6406250, 0.5526316, 0.5454545, NA, 0.0000000)
Slc4a1_naive_FDR = c(NA, 0.69333333, 0.03773585, 0.00000000, NA, 0.31067961)

Cd4_adj_FDR = c(NA, 0.3500000, 0.1046512, 0.0000000, NA, 0.7575758)
Rps19_adj_FDR = c(NA, 0.5806452, 0.3362069, 0.3238095, NA, 0.8987342)
Slc4a1_adj_FDR = c(NA, 0.7111111, 0.5497076, 0.4033613, NA, 0.4285714)

# get FDR mat
FDR_mat_naive = cbind(Cd4_naive_FDR, Rps19_naive_FDR, Slc4a1_naive_FDR)
FDR_mat_adj = cbind(Cd4_adj_FDR, Rps19_adj_FDR, Slc4a1_adj_FDR)
colnames(FDR_mat_naive) = c('Cd4', 'Rps19', 'Slc4a1')
colnames(FDR_mat_adj) = c('Cd4', 'Rps19', 'Slc4a1')

# plot FDR bar plot
for (NMF_F_n in 1:dim(FDR_mat_naive)[1]){
pdf(paste0('NMF_F',NMF_F_n,'.FDR.pdf'), width = 4, height = 8)
par(mfrow=c(2,1))
barplot(FDR_mat_naive[NMF_F_n,], ylim=c(0,1), ylab='FDR', xlab='NMF factor', main='naive')
box()
barplot(FDR_mat_adj[NMF_F_n,], ylim=c(0,1), ylab='FDR', xlab='NMF factor', main='adjusted')
box()
dev.off()
}





FDR_p_thresh = get_FDR_p_thresh(nmf_result_0_W[,j], nmf_result_W_mat_long[,j], 0.05, 1e-10, -0.0001)

0.0126 0.0126 0.0016 0.0020 0.0126 0.0126
0.0068 0.0185 0.0126 0.0203 0.0126 0.0126

Cd4
NaN 0.50000000 0.02666667 0.00000000        NaN 0.32380952
NaN 0.4444444 0.0700000 0.0000000       NaN 1.0000000
NaN 0.53333333 0.04255319 0.01315789        NaN 0.47368421

Rps19
NaN 0.7500000 0.2962963 0.2183908       NaN 0.6802326
NaN 0.7368421 0.2954545 0.2589286       NaN 1.0000000

Slc4a1
NaN 0.8571429 0.5539906 0.1168831       NaN 0.2028986
NaN 0.8529412 0.5485437 0.4275862       NaN       NaN



reconstruct_mat_0 = (nmf_result_0_W*nmf_result_0_W_binary) %*% nmf_result_0_H
reconstruct_mat_bg = (nmf_result_W*nmf_result_W_binary) %*% nmf_result_H

# plot NMFs_mm10 binary heatmap
my_colorbar_binary=colorRampPalette(c('white','black'))(n = 100)
png(paste0(output_folder, '/', hg38_gene, '.', mm10_gene, '.cor.heatmap.png.cor.mat.txt.NMFs_mm10.reconstruct_heatmap.png'), width = 5000, height = 5000)
pheatmap(reconstruct_mat_0, cluster_rows=F, cluster_cols=F, width = 2, height = 6, legend=F, show_rownames=F, show_colnames=F)
dev.off()

my_colorbar_binary=colorRampPalette(c('white','black'))(n = 100)
png(paste0(output_folder, '/', hg38_gene, '.', mm10_gene_bg, '.cor.heatmap.png.cor.mat.txt.NMFs_mm10.reconstruct_heatmap.png'), width = 5000, height = 5000)
pheatmap(reconstruct_mat_bg, cluster_rows=F, cluster_cols=F, width = 2, height = 6, legend=F, show_rownames=F, show_colnames=F)
dev.off()










# cluster factor matrix
nmf_W_H = rbind(nmf_result_W, t(nmf_result_H))
nmf_W_H_order_tree = hclust(dist(t(nmf_W_H)), method = 'average')
nmf_W_H_order = nmf_W_H_order_tree$order


# plot NMF factor i bar plot
for (i in 1:dim(nmf_result_W)[2]){
    # plot hg38
    png(paste0(output_folder, '/', hg38_gene, '.', mm10_gene_bg, '.cor.heatmap.png.cor.mat.txt.NMF.hg38.', i,'.bar.png'), width = 1000, height = 300)
    cor_factor = nmf_result_H[i,]
    span = 0.03
    cor_factor_sm = loess(cor_factor ~ c(1:dim(nmf_result_H)[2]), span=span)$fitted
    plot(1:dim(nmf_result_H)[2], cor_factor_sm, type='h', ylim=c(0,1))
    dev.off()
    # plot mm10
    png(paste0(output_folder, '/', hg38_gene, '.', mm10_gene_bg, '.cor.heatmap.png.cor.mat.txt.NMF.mm10.', i,'.bar.png'), width = 1000, height = 300)
    cor_factor = nmf_result_W[,i]
    cor_factor_sm = loess(cor_factor ~ c(1:dim(nmf_result_W)[1]), span=span)$fitted
    plot(1:dim(nmf_result_W)[1], cor_factor_sm, type='h', ylim=c(0,1))
    dev.off()
}


# plot NMF factor i z score p-value bar plot
for (i in 2:dim(nmf_result_W)[2]){
    # plot hg38
    png(paste0(output_folder, '/', hg38_gene, '.', mm10_gene_bg, '.cor.heatmap.png.cor.mat.txt.NMF.hg38.', i,'.zp.bar.png'), width = 1000, height = 100)
    par(mar = rep(0.5, 4))
    cor_factor = nmf_result_H[i,]
    span = 0.03
    cor_factor_sm = loess(cor_factor ~ c(1:dim(nmf_result_H)[2]), span=span)$fitted
    cor_factor_sm[cor_factor_sm<=0] = 0.1
    cor_factor_sm_zp = -log10(zp2(cor_factor_sm, 0.01))
    #plot(1:dim(nmf_result_H)[2], cor_factor_sm_zp, type='h', ylim=c(0,16), xlim=c(1,dim(nmf_result_H)[2]), xlab='', ylab='', xaxt='n', yaxt='n')
    barplot(cor_factor_sm_zp, xlab='', ylab='', xaxt='n', yaxt='n', ylim=c(0,8))
    dev.off()
    # fg
    cor_factor = nmf_result_0_W[,i]
    cor_factor_sm = loess(cor_factor ~ c(1:dim(nmf_result_0_W)[1]), span=span)$fitted
    cor_factor_sm[cor_factor_sm<=0] = 1e-16
    cor_factor_bg = nmf_result_W[,i]
    cor_factor_sm_bg = loess(cor_factor_bg ~ c(1:dim(nmf_result_W)[1]), span=span)$fitted
    cor_factor_sm_bg[cor_factor_sm_bg<=0] = 1e-16
    cor_factor_sm_zp_raw = zp2_bg(cor_factor_sm, cor_factor_sm_bg)
    cor_factor_sm_zp_raw_FDR = p.adjust(cor_factor_sm_zp_raw, method = 'fdr')
    print(sum(cor_factor_sm_zp_raw_FDR<0.1))
    # plot mm10 bg
    png(paste0(output_folder, '/', hg38_gene, '.', mm10_gene_bg, '.cor.heatmap.png.cor.mat.txt.NMF.mm10.', i,'.zp.bar.png'), width = 1000, height = 100)
    par(mar = rep(0.5, 4))
    cor_factor = nmf_result_W[,i]
    cor_factor_sm = loess(cor_factor ~ c(1:dim(nmf_result_W)[1]), span=span)$fitted
    cor_factor_sm[cor_factor_sm<=0] = 1e-16
    cor_factor_sm_zp_raw = zp2(cor_factor_sm, 0.01)
    cor_factor_sm_zp_raw_FDR = p.adjust(cor_factor_sm_zp_raw, method = 'fdr')
    print(sum(cor_factor_sm_zp_raw_FDR<0.1))
    cor_factor_sm_zp = -log10(cor_factor_sm_zp_raw)
    cor_factor_sm_zp[!is.finite(cor_factor_sm_zp)] = 16
    # plot barplot without axis labels
    #plot(1:dim(nmf_result_W)[1], cor_factor_sm_zp, type='h', ylim=c(0,16), xlim=c(1,dim(nmf_result_W)[1]), xlab='', ylab='', xaxt='n', yaxt='n')
    barplot(cor_factor_sm_zp, xlab='', ylab='', xaxt='n', yaxt='n', ylim=c(0,8))
    dev.off()
}

# plot NMF factor i z score p-value bar plot bg gene adjusted
for (i in 2:dim(nmf_result_W)[2]){
    # mm10
    # fg
    cor_factor = nmf_result_0_W[,i]
    cor_factor_sm = loess(cor_factor ~ c(1:dim(nmf_result_0_W)[1]), span=span)$fitted
    cor_factor_sm[cor_factor_sm<=0] = 1e-16
    # bg
    cor_factor_bg = nmf_result_W[,i]
    cor_factor_sm_bg = loess(cor_factor_bg ~ c(1:dim(nmf_result_W)[1]), span=span)$fitted
    cor_factor_sm_bg[cor_factor_sm_bg<=0] = 1e-16
    # zp with bg gene adjusted
    cor_factor_sm_zp_raw = zp2_bg(cor_factor_sm, cor_factor_sm_bg)
    cor_factor_sm_zp_raw_FDR = p.adjust(cor_factor_sm_zp_raw, method = 'fdr')
    print(sum(cor_factor_sm_zp_raw_FDR<0.1))
}

Gata1
[1] 23
[1] 51
[1] 35
[1] 0
[1] 71

Cd4
[1] 34
[1] 29
[1] 50
[1] 0
[1] 36

Rps19
[1] 41
[1] 63
[1] 42
[1] 0
[1] 0

Slc4a1
[1] 52
[1] 2
[1] 0
[1] 0
[1] 32

[1] 3
[1] 43
[1] 26
[1] 0
[1] 70


apply(nmf_result_W_mat_long,2,function(x) mean(x[x!=0]))
[1] 0.15533594 0.07678033 0.03423006 0.07439433 0.32601713 0.04284531

> apply(nmf_result_W_mat_long,2,function(x) mean(x[x!=0]))
[1] 0.15413623 0.07678739 0.03397224 0.06475345 0.31497731 0.08751369

apply(nmf_result_W_mat_long,2,function(x) sum(x!=0))
[1] 1102 1102 1101 1099 1102 1082


Cd4
NaN 0.5964912 0.3625000 0.5882353       NaN 0.3364486
Rps19
NaN 0.6406250 0.5526316 0.5454545       NaN 0.0000000
Slc4a1
NaN 0.69333333 0.03773585 0.00000000        NaN 0.31067961

Cd4
NaN 0.4444444 0.0700000 0.0000000       NaN 1.0000000
Rps19
NaN 0.7368421 0.2954545 0.2589286       NaN 1.0000000
Slc4a1
NaN 0.8529412 0.5485437 0.4275862       NaN       NaN





p.adjust = function (p, method = p.adjust.methods, n = length(p)) 
{
    method <- match.arg(method)
    if (method == "fdr") 
        method <- "BH"
    nm <- names(p)
    p <- as.numeric(p)
    p0 <- setNames(p, nm)
    if (all(nna <- !is.na(p))) 
        nna <- TRUE
    else p <- p[nna]
    lp <- length(p)
    stopifnot(n >= lp)
    if (n <= 1) 
        return(p0)
    if (n == 2 && method == "hommel") 
        method <- "hochberg"
    p0[nna] <- switch(method, bonferroni = pmin(1, n * p), holm = {
        i <- seq_len(lp)
        o <- order(p)
        ro <- order(o)
        pmin(1, cummax((n + 1L - i) * p[o]))[ro]
    }, hommel = {
        if (n > lp) p <- c(p, rep.int(1, n - lp))
        i <- seq_len(n)
        o <- order(p)
        p <- p[o]
        ro <- order(o)
        q <- pa <- rep.int(min(n * p/i), n)
        for (j in (n - 1L):2L) {
            ij <- seq_len(n - j + 1L)
            i2 <- (n - j + 2L):n
            q1 <- min(j * p[i2]/(2L:j))
            q[ij] <- pmin(j * p[ij], q1)
            q[i2] <- q[n - j + 1L]
            pa <- pmax(pa, q)
        }
        pmax(pa, p)[if (lp < n) ro[1L:lp] else ro]
    }, hochberg = {
        i <- lp:1L
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin((n + 1L - i) * p[o]))[ro]
    }, BH = {
        i <- lp:1L
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin(n/i * p[o]))[ro]
    }, BY = {
        i <- lp:1L
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        q <- sum(1/(1L:n))
        pmin(1, cummin(q * n/i * p[o]))[ro]
    }, none = p)
    p0
}