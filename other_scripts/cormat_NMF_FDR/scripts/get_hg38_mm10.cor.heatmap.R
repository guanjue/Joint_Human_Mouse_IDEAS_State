#library(DescTools)
args = commandArgs(trailingOnly=TRUE)

hg38_gene_state = args[1]
mm10_gene_state = args[2]
output_file = args[3]
state_file = args[4]

#hg38_gene_state = 'hg38.gene.GATA1.matched_ct.state.bed'
#mm10_gene_state = 'mm10.gene.Gata1.matched_ct.state.bed'
#state_file = 'IDEAS.EpigeneticState.mean_signal_mat.txt'
#output_file = 'GATA1.cor.heatmap.png'


### read state signal mat
state_mat_od = as.matrix(read.table(state_file, header=T, sep='\t'))

### state gene assignment
d1 = read.table(hg38_gene_state, header=F, sep='\t', comment.char='~')
d2 = read.table(mm10_gene_state, header=F, sep='\t', comment.char='~')

### get state
d1s = d1[,-c(1:3)]
d2s = d2[,-c(1:3)]

### state to signal matrix
d1s_sigmat = matrix(0, nrow=dim(d1s)[1], ncol=dim(d1s)[2]*dim(state_mat_od)[2])
d2s_sigmat = matrix(0, nrow=dim(d2s)[1], ncol=dim(d2s)[2]*dim(state_mat_od)[2])
###
for (i in 1:dim(d1s)[1]){
	sig_i = c(state_mat_od[unlist(d1s[i,]+1),])
	d1s_sigmat[i,] = as.numeric(state_mat_od[unlist(d1s[i,]+1),])
}
###
for (i in 1:dim(d2s)[1]){
	sig_i = c(state_mat_od[unlist(d2s[i,]+1),])
	d2s_sigmat[i,] = as.numeric(state_mat_od[unlist(d2s[i,]+1),])
}
###
rownames(d1s_sigmat) = paste('H',1:dim(d1s_sigmat)[1], sep=':')
rownames(d2s_sigmat) = paste('M',1:dim(d2s_sigmat)[1], sep=':')

###
d1s_sigmat_rowMax = log(apply(d1s_sigmat,1,max)+1)
d1s_sigmat_rowMax = d1s_sigmat_rowMax/max(d1s_sigmat_rowMax)
d2s_sigmat_rowMax = log(apply(d2s_sigmat,1,max)+1)
d2s_sigmat_rowMax = d2s_sigmat_rowMax/max(d2s_sigmat_rowMax)

### get cor matrix
d12_cor_mat = as.matrix(cor(t(d2s_sigmat+matrix(rnorm(dim(d2s_sigmat)[1]*dim(d2s_sigmat)[2], sd=0.2), nrow=dim(d2s_sigmat)[1], ncol=dim(d2s_sigmat)[2]) ), t(d1s_sigmat+matrix(rnorm(dim(d1s_sigmat)[1]*dim(d1s_sigmat)[2], sd=0.2), nrow=dim(d1s_sigmat)[1], ncol=dim(d1s_sigmat)[2]) )))
d12_cor_mat[is.na(d12_cor_mat)] = 0

d12_cor_mat_adj = t(apply(d12_cor_mat,1,function(x) x*1))
d12_cor_mat_adj = apply(d12_cor_mat_adj,2,function(x) x*1)

# write out cor matrix
write.table(d12_cor_mat_adj, paste(output_file, '.cor.mat.txt', sep=''), row.names = F, col.names = F, quote = F, sep = '\t')

library(pheatmap)
breaksList = seq(-1, 1, by = 0.001)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList))
png(output_file, height=800, width=800)
colnames(d12_cor_mat_adj) = NULL
rownames(d12_cor_mat_adj) = NULL
pheatmap(d12_cor_mat_adj, cluster_rows=F, cluster_cols=F, color=my_colorbar, breaks = breaksList)
dev.off()
