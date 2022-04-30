cts = c("B", "CLP", "CMP", "EOS", "ERY", "GMP", "LSK", "MK", "MONc", "MONp", "MPP", "NEU", "NK", "CD4", "CD8")
cts_mouse = c("CFUE", "CFUMK", "CMP", "ERYfl", "GMP", "LSK", "MEP", "MONO", "NEU", "iMK", "G1E", "ER4")


### human
sig0 = read.table(paste('coe_score_no',cts[1],'/statep_rna_coe_heatmap.human.all.ccre.withcorfilter.txt', sep=''), header=T)

for (i in 2:length(cts)){
	sig0_i = read.table(paste('coe_score_no',cts[i],'/statep_rna_coe_heatmap.human.all.ccre.withcorfilter.txt', sep=''), header=T)
	sig0 = sig0+sig0_i
}

sig0 = sig0/length(cts)


### mouse
sig0_mouse = read.table(paste('../coe_analysis_mouse/coe_score_no',cts_mouse[1],'/statep_rna_coe_heatmap.mouse.all.ccre.withcorfilter.txt', sep=''), header=T)

for (i in 2:length(cts_mouse)){
	sig0_i = read.table(paste('../coe_analysis_mouse/coe_score_no',cts_mouse[i],'/statep_rna_coe_heatmap.mouse.all.ccre.withcorfilter.txt', sep=''), header=T)
	sig0_mouse = sig0_mouse+sig0_i
}

sig0_mouse = sig0_mouse/length(cts_mouse)

sig0_mouse_min_non0 = min(sig0_mouse[sig0_mouse>0])
sig0_mouse = sig0_mouse-sig0_mouse_min_non0
sig0_mouse = sig0_mouse /max(sig0_mouse) * (max(sig0)-sig0_mouse_min_non0)
sig0_mouse = sig0_mouse+sig0_mouse_min_non0




library(pheatmap)


sig0_ave = (sig0+sig0_mouse)/2
#sig0_ave_direction = ((sig0+sig0_mouse)>=0)*1
#print(sig0_ave_direction)
#sig0_ave_direction[sig0_ave_direction==0] = -1
#print(sig0_ave_direction)
#sig0_ave = sqrt(abs(sig0*sig0_mouse)) * sig0_ave_direction

plot_lim_P = max(abs(sig0_ave[,1]))
plot_lim_D = max(abs(sig0_ave[,2]))
plot_lim_PD = max(abs(sig0_ave))
pdf(paste('statep_rna_coe_heatmap.HM.all.ccre.withcorfilter.AVE.pdf', sep=''), width=3)
rank = c(2,1,4,3,6,10,11,5,9,13,8,19,12,25,14,23,22,20,21,17,18,7,16,15,24)
breaksList = seq(-plot_lim_PD, plot_lim_PD, by = 0.001)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList))
col_breaks = c(seq(0, 2000,length=33))
pheatmap(sig0_ave[rank,], color=my_colorbar, breaks = breaksList, cluster_cols = FALSE,cluster_rows=F, clustering_method = 'average',annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()

pdf(paste('statep_rna_coe_heatmap.HM.P.all.ccre.withcorfilter.AVE.pdf', sep=''), width=1)
rank = c(2,1,4,3,6,10,11,5,9,13,8,19,12,25,14,23,22,20,21,17,18,7,16,15,24)
breaksList = seq(-plot_lim_P, plot_lim_P, by = 0.001)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList))
col_breaks = c(seq(0, 2000,length=33))
pheatmap(cbind(sig0_ave[rank,1]), color=my_colorbar, breaks = breaksList, cluster_cols = FALSE,cluster_rows=F, clustering_method = 'average',annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()

pdf(paste('statep_rna_coe_heatmap.HM.D.all.ccre.withcorfilter.AVE.pdf', sep=''), width=1)
rank = c(2,1,4,3,6,10,11,5,9,13,8,19,12,25,14,23,22,20,21,17,18,7,16,15,24)
breaksList = seq(-plot_lim_D, plot_lim_D, by = 0.001)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList))
col_breaks = c(seq(0, 2000,length=33))
pheatmap(cbind(sig0_ave[rank,2]), color=my_colorbar, breaks = breaksList, cluster_cols = FALSE,cluster_rows=F, clustering_method = 'average',annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()

pdf(paste('statep_rna_coe_heatmap.HM.PDsum.all.ccre.withcorfilter.AVE.pdf', sep=''), width=1)
rank = c(2,1,4,3,6,10,11,5,9,13,8,19,12,25,14,23,22,20,21,17,18,7,16,15,24)
breaksList = seq(-plot_lim_PD, plot_lim_PD, by = 0.001)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList))
col_breaks = c(seq(0, 2000,length=33))
pheatmap(cbind(sig0_ave[rank,1]+sig0_ave[rank,2]), color=my_colorbar, breaks = breaksList, cluster_cols = FALSE,cluster_rows=F, clustering_method = 'average',annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()

write.table(sig0_ave, 'statep_rna_coe_heatmap.HM.all.ccre.withcorfilter.AVE.txt', quote=F, col.names=T, row.names=T, sep='\t')


pdf('coe_mouse_vs_human.pdf', width=3.6, height=4)
plot_lim_max = max(c(max(abs(sig0)), max(abs(sig0_mouse))) )
plot_lim_min = min(c(min((sig0)), min((sig0_mouse))) )
plot(as.numeric(as.matrix(sig0)), as.numeric(as.matrix(sig0_mouse)), xlim=c(plot_lim_min,plot_lim_max), ylim=c(plot_lim_min,plot_lim_max), xlab='Human', ylab='Mouse', 
	main = paste('Pearson.Cor = ', round(cor(as.numeric(as.matrix(sig0)), as.numeric(as.matrix(sig0_mouse))), 3), sep=''))
abline(0,1, col='red')
dev.off()


