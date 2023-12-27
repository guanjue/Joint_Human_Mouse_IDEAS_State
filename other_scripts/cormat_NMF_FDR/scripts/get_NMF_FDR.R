args = commandArgs(trailingOnly=TRUE)

working_dir = args[1]
hg38_gene = args[2]
mm10_gene = args[3]
random_picked_bg_genes_file_mm10 = args[4]
random_picked_bg_genes_file_hg38 = args[5]
fdr_thresh = as.numeric(args[6])
n_comp = as.numeric(args[7])

library(pheatmap)

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

# p-vlaue of Z score based on a background distribution 
zp2_bg = function(x, xbg) {
    # 1st round z score
    z = (x - mean(xbg)) / sd(xbg)
    # get z score one-side p value
    zp = 1 - pnorm(z)
    return(zp)
}

# get p-vlaue of Z score based on a background distribution with smoothing
zp2_sm_bg = function(x, xbg, span) {
    # read x
    x_sm = loess(x ~ c(1:length(x)), span=span)$fitted
    x_sm[x_sm<=0] = min(x_sm[x_sm>0])
    x_sm_zp = zp2_bg(x_sm, xbg)
    return(x_sm_zp)
}

# get p-value of Z score equivalent to FDR = fdr_thresh
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
    # check if FDR_p_thresh is defined
    if (!exists('FDR_p_thresh')) {
        FDR_p_thresh = 0
    }
    return(FDR_p_thresh)
}

# get background genes NMF matrix function
background_gene_NMF_mat = function(random_picked_bg_genes_file, input_folder_cormat, target_gene, bed_mm10, bed_hg38, nmf_result_0, transpose = FALSE){
    # read random_pick_bg_genes
    random_picked_bg_genes = read.table(random_picked_bg_genes_file, header = F)

    # read cor_mat for random picked bg genes
    nmf_result_W_mat_long = c()

    for (i in 1:nrow(random_picked_bg_genes)) {
        # read gene name
        mm10_gene_bg_i = random_picked_bg_genes[i,1]
        print(mm10_gene_bg_i)
        # get mm10_gene_bg file
        input_correlation_mat_bg_i = paste0(input_folder_cormat, target_gene, '.', mm10_gene_bg_i, '.cor.heatmap.png.cor.mat.txt')

        if (file.exists(input_correlation_mat_bg_i)) {
        # read correlation matrix
        d_mk_cor_bg_i = as.matrix(read.table(input_correlation_mat_bg_i, header=F, sep='\t', comment.char='~'))
        d_mk_cor_mat_bg_i = as.matrix(d_mk_cor_bg_i)

        # check if the dimension of d_mk_cor_mat_bg_i is correct
        if ((dim(d_mk_cor_mat_bg_i)[1] != dim(bed_mm10)[1]) | (dim(d_mk_cor_mat_bg_i)[2] != dim(bed_hg38)[1])) {
          if (((dim(d_mk_cor_mat_bg_i)[1] - dim(bed_mm10)[1]) == 1) & (dim(d_mk_cor_mat_bg_i)[2] == dim(bed_hg38)[1])) {
            d_mk_cor_mat_bg_i = d_mk_cor_mat_bg_i[-1,]
          } else if (((dim(d_mk_cor_mat_bg_i)[2] - dim(bed_hg38)[1]) == 1) & (dim(d_mk_cor_mat_bg_i)[1] == dim(bed_mm10)[1])) {
            d_mk_cor_mat_bg_i = d_mk_cor_mat_bg_i[,-1]
          }
            else if (((dim(d_mk_cor_mat_bg_i)[1] - dim(bed_mm10)[1]) == 1) & ((dim(d_mk_cor_mat_bg_i)[2] - dim(bed_hg38)[1]) == 1)) {
                d_mk_cor_mat_bg_i = d_mk_cor_mat_bg_i[-1,-1]
          } else {
            next
          }
        }
        # add rownames and colnames
        colnames(d_mk_cor_mat_bg_i) = apply(bed_hg38, 1, function(x) paste0('H_', x[1], '_', x[2], '_', x[3]))
        rownames(d_mk_cor_mat_bg_i) = apply(bed_mm10, 1, function(x) paste0('M_', x[1], '_', x[2], '_', x[3]))

        # get positive correlation matrix
        d_mk_cor_mat_bg_pos_i = d_mk_cor_mat_bg_i
        d_mk_cor_mat_bg_pos_i[d_mk_cor_mat_bg_pos_i<0] = 0

        # nmf fix H
        set.seed(2019)
        if (transpose == FALSE) {
            nmf_result_i <- nmf_fix_H(d_mk_cor_mat_bg_pos_i, H = nmf_result_0)
        } else {
            nmf_result_i <- nmf_fix_H(d_mk_cor_mat_bg_pos_i, H = t(nmf_result_0))
        }

        # Extract the independent components (ICs) and the mixing matrix
        nmf_result_W_i <- nmf_result_i$W

        # add to nmf_result_W_mat_long
        nmf_result_W_mat_long = rbind(nmf_result_W_mat_long, nmf_result_W_i)
    }}
    return(nmf_result_W_mat_long)
}

# get binary matrix based on random selected bg genes adjusted FDR for foreground gene
get_bggene_adj_FDR_binary_mat = function(nmf_result_0_W, nmf_result_W_mat_long, fdr_thresh){
    nmf_result_0_W_binary = nmf_result_0_W
    FDR_p_thresh_0_vec = c()
    for (j in 1:dim(nmf_result_0_W)[2]){
        print(j)
        # get nmf_result_0_W_j and nmf_result_W_mat_long_j
        nmf_result_0_W_j = nmf_result_0_W[,j]
        nmf_result_W_mat_long_j = nmf_result_W_mat_long[,j]
        # get FDR_p_thresh
        FDR_p_thresh_0_j = get_FDR_p_thresh(nmf_result_0_W_j, nmf_result_W_mat_long_j, 0.05, 1e-10, -0.0001, fdr_thresh)
        FDR_p_thresh_0_vec = c(FDR_p_thresh_0_vec, FDR_p_thresh_0_j)
        # get binary matrix based on p-value equivalent to FDR = fdr_thresh
        nmf_result_0_W_sm_zpfdr = (zp2_sm_bg(nmf_result_0_W_j, nmf_result_W_mat_long_j, 0.03) < FDR_p_thresh_0_j) *1
        nmf_result_0_W_binary[,j] = nmf_result_0_W_sm_zpfdr
    }
    return(nmf_result_0_W_binary)
}

# plot NMFs binary heatmap
plot_NMF_FDR_binary_heatmap = function(nmf_result_0_W_binary, output_folder, outputname){
    my_colorbar_binary=colorRampPalette(c('white','black'))(n = 100)
    png(paste0(output_folder, '/', outputname), width = 300, height = 1000)
    pheatmap(nmf_result_0_W_binary, cluster_rows=F, cluster_cols=F, width = 2, height = 6, legend=F, show_rownames=F, show_colnames=F, color=my_colorbar_binary)
    dev.off()
}

##############################################################################################################

# set working dir
setwd(working_dir)

# output plot folder
output_folder = paste0('NMF_reconstruction_', hg38_gene, '_', mm10_gene)
input_correlation_mat = paste0(hg38_gene, '.', mm10_gene, '.cor.heatmap.png.cor.mat.txt')

# read coordinates
hg38_coord = paste0('hg38.gene.', hg38_gene, '.matched_ct.state.bed')
mm10_coord = paste0('mm10.gene.', mm10_gene, '.matched_ct.state.bed')

# read correlation matrix
d_mk_cor = as.matrix(read.table(input_correlation_mat, header=F, sep='\t', comment.char='~'))
d_mk_cor_mat = as.matrix(d_mk_cor)

# read coordinates
bed_hg38 = read.table(hg38_coord, header=F, sep='\t', comment.char='~')
bed_mm10 = read.table(mm10_coord, header=F, sep='\t', comment.char='~')
# add rownames and colnames
colnames(d_mk_cor_mat) = apply(bed_hg38, 1, function(x) paste0('H_', x[1], '_', x[2], '_', x[3]))
rownames(d_mk_cor_mat) = apply(bed_mm10, 1, function(x) paste0('M_', x[1], '_', x[2], '_', x[3]))

# get positive correlation matrix
d_mk_cor_mat_pos = d_mk_cor_mat
d_mk_cor_mat_pos[d_mk_cor_mat_pos<0] = 0

# mkdir
dir.create(output_folder)

# decide the number of components based on the BIC
# check if the NMFs.num.BIC.png file not exists
if (!file.exists(paste0(output_folder, '/', hg38_gene, '.', mm10_gene, '.NMFs.num.BIC.pdf'))) {
    # decide the number of components based on the BIC
    set.seed(2023)
    BIC_mat = c()
    NMF_n_max = 20
    for (j in 1:1){
        print(j)
        BIC_vec = c()
        for (i in seq(1, NMF_n_max, by = 1)) {
            print(i)
            nmf_result_tmp <- nmf(d_mk_cor_mat_pos, k = i)
            W <- nmf_result_tmp$W
            H <- nmf_result_tmp$H
            reconstructed_mat <- W %*% H
            # get BIC of the model
            BIC = nrow(d_mk_cor_mat_pos) * ncol(d_mk_cor_mat_pos) * log(sum((reconstructed_mat - d_mk_cor_mat_pos)^2)) + i * (nrow(d_mk_cor_mat_pos) + ncol(d_mk_cor_mat_pos)) * log(nrow(d_mk_cor_mat_pos) * ncol(d_mk_cor_mat_pos))
            print(BIC)
            BIC_vec = c(BIC_vec, BIC)
            # save i = 20 reconstructed_mat
            if (i == NMF_n_max) {
                reconstructed_mat_W = W %*% t(W)
                reconstructed_mat_H = t(H) %*% H
            }
        }
        BIC_mat = rbind(BIC_mat, BIC_vec)
    }
    # plot BIC
    pdf(paste0(output_folder, '/', hg38_gene, '.', mm10_gene, '.NMFs.num.BIC.pdf'))
    plot(seq(1, NMF_n_max, by = 1), BIC_mat[1,], type='l', xlab='Number of components', ylab='BIC', main='BIC')
    points(seq(1, NMF_n_max, by = 1), BIC_mat[1,])
    dev.off()
}

#  NMF decomposition
set.seed(2019)
#n_comp = 6
nmf_result_0 <- nmf(d_mk_cor_mat_pos, k = n_comp)

# Extract the independent components (ICs) and the mixing matrix
nmf_result_0_W <- nmf_result_0$W
nmf_result_0_H <- nmf_result_0$H

# cluster factor matrix
nmf_W_H = rbind(nmf_result_0_W, t(nmf_result_0_H))
nmf_W_H_order_tree = hclust(dist(t(nmf_W_H)), method = 'average')
nmf_W_H_order = nmf_W_H_order_tree$order

# plot NMFs_mm10 heatmap
png(paste0(output_folder, '/', mm10_gene, '.cor.heatmap.png.cor.mat.txt.NMFs.mm10.heatmap.png'), width = 300, height = 1000)
pheatmap(nmf_result_0_W, cluster_rows=F, cluster_cols=F, width = 2, height = 6, legend=F, show_rownames=F, show_colnames=F)
dev.off()
# plot NMFs_hg38 heatmap
png(paste0(output_folder, '/', hg38_gene, '.cor.heatmap.png.cor.mat.txt.NMFs.hg38.heatmap.png'), width = 300, height = 1000)
pheatmap(t(nmf_result_0_H), cluster_rows=F, cluster_cols=F, width = 2, height = 6, legend=F, show_rownames=F, show_colnames=F)
dev.off()

##############################################################################################################

# get background genes NMF matrix mm10
nmf_result_W_mat_long_mm10 = background_gene_NMF_mat(random_picked_bg_genes_file_mm10, 'randomly_pick_bg_genes_mm10/', hg38_gene, bed_mm10, bed_hg38, nmf_result_0_H)
print(dim(nmf_result_W_mat_long_mm10))
print(dim(nmf_result_0_W))
# get NMF FDR binary matrix mm10
nmf_result_0_W_binary_mm10 = get_bggene_adj_FDR_binary_mat(nmf_result_0_W, nmf_result_W_mat_long_mm10, fdr_thresh)
# write binary matrix
write.table(nmf_result_0_W_binary_mm10, paste0(output_folder, '/', mm10_gene, '.cor.mat.txt.NMFs.mm10.binary_FDRbgadj.txt'), sep='\t', quote=F, row.names=T, col.names=F)
# plot NMFs_mm10 binary heatmap
plot_NMF_FDR_binary_heatmap(nmf_result_0_W_binary_mm10, output_folder, paste0(mm10_gene, '.cor.mat.txt.NMFs.mm10.binary_FDRbgadj.heatmap.png'))

# get background genes NMF matrix hg38
nmf_result_W_mat_long_hg38 = background_gene_NMF_mat(random_picked_bg_genes_file_hg38, 'randomly_pick_bg_genes_hg38/', mm10_gene, bed_hg38, bed_mm10, nmf_result_0_W, transpose = TRUE)
print(dim(nmf_result_W_mat_long_hg38))
print(dim(nmf_result_0_H))
# get NMF FDR binary matrix hg38
nmf_result_0_W_binary_hg38 = get_bggene_adj_FDR_binary_mat(t(nmf_result_0_H), nmf_result_W_mat_long_hg38, fdr_thresh)
# write binary matrix
write.table(nmf_result_0_W_binary_hg38, paste0(output_folder, '/', hg38_gene, '.cor.mat.txt.NMFs.hg38.binary_FDRbgadj.txt'), sep='\t', quote=F, row.names=T, col.names=F)
# plot NMFs_hg38 binary heatmap
plot_NMF_FDR_binary_heatmap(nmf_result_0_W_binary_hg38, output_folder, paste0(hg38_gene, '.cor.mat.txt.NMFs.hg38.binary_FDRbgadj.heatmap.png'))

##############################################################################################################

# plot cross factor correlation matrix
# get top quantile of NMFs
nmf_result_0_W_binary_top_quantile = apply(nmf_result_0_W, 2, function(x) x > quantile(x, 0.9))*1
nmf_result_0_H_binary_top_quantile = apply(t(nmf_result_0_H), 2, function(x) x > quantile(x, 0.9))*1

# get mean correlation between different NMF factors
cross_factor_cor_mat = matrix(0, nrow = dim(nmf_result_0_W_binary_hg38)[2], ncol = dim(nmf_result_0_W_binary_mm10)[2])
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
colnames(cross_factor_cor_mat) = paste0(hg38_gene, '_', 1:dim(nmf_result_0_W_binary_hg38)[2])
rownames(cross_factor_cor_mat) = paste0(mm10_gene, '_', 1:dim(nmf_result_0_W_binary_mm10)[2])
# my colorbar white to red
my_colorbar = colorRampPalette(c('white','red'))(n = 100)
pdf(paste0(output_folder, '/', hg38_gene, '.', mm10_gene, '.cross_factor_cor_mat.pdf'), width = 10, height = 10)
pheatmap(cross_factor_cor_mat, cluster_rows=F, cluster_cols=F, color = my_colorbar, legend=T, show_rownames=T, show_colnames=T, cex=1.5)
dev.off()

# write cross_factor_cor_mat
write.table(cross_factor_cor_mat, paste0(output_folder, '/', hg38_gene, '.', mm10_gene, '.cross_factor_cor_mat.txt'), sep='\t', quote=F, row.names=T, col.names=T)




