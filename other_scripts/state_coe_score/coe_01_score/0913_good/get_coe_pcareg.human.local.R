R
library(data.table)
library(LSD)
library(pheatmap)
library(pastecs)
library(lmvar)
library(mclust)

R2 = function(obs, pred){
	r2 = 1-mean((obs-pred)^2)/mean((obs-mean(obs))^2)
	return(r2)
}

R2_0 = function(obs, pred){
	r2 = 1-mean((obs-pred)^2)/mean((obs-0)^2)
	return(r2)
}

quantile_norm = function(x){
	xm = rowMeans(x)
	refsig_sort = xm[order(xm)]
	for (i in 1:dim(x)[2]){
		sigtmp = x[,i]
		x[,i] = refsig_sort[rank(sigtmp)]
	}
	return(x)
}

setwd('/gpfs/group/yzz2/default/legacy/group/projects/vision_joint_human_mouse/coefficients_human')

sp_list = c('B', 'B', 'CD34', 'CD34', 'CD4', 'CD4', 'CD8', 'CD8', 'CLP', 'CLP', 'CMP', 'CMP', 'EOS', 'ERY', 'ERY', 'GMP', 'GMP', 'H1', 'H1', 'H2', 'H2', 'LMPP', 'LMPP', 'LSK', 'LSK', 'MK', 'MK', 'MONc', 'MONc', 'MONp', 'MONp', 'MPP', 'MPP', 'NEU', 'NEU', 'NK', 'NK', 'PBMC', 'PBMC')
rna_list = c('LSK', 'LSK', 'ERY', 'ERY', 'CD4', 'CD4', 'CD8', 'CD8', 'B', 'B', 'CMP', 'CMP', 'MONp', 'MONp', 'NEU', 'NEU', 'MONc', 'MONc', 'GMP', 'GMP', 'CFUE', 'NK', 'NK', 'MK', 'MK', 'CLP', 'CLP', 'MPP', 'MPP', 'EOS', 'EOS')
common_ct = intersect(sp_list, rna_list)

smallnum = 1e-1
smallnum_sp = 1
plot_lim_all = c(-2.5,5)
set.seed(2019)

### get RNA
d = read.table('HumanVISION_RNAseq_hg38_genes_tpm.idsort.protein_coding.txt', header=F)
ds = d[,-c(1:6)]

### QTnorm
ds_qt = quantile_norm(ds)

### split gene groups
dms = rowMeans(ds_qt)
dms_log = log(dms+smallnum)
mod = densityMclust(dms_log,G=2)
gmm_thresh = max(dms_log[mod$classification==1])
rna_l = mod$classification==1
rna_h = mod$classification==2
pdf('dms_log_rna.hist.pdf', width=5, height=5)
plot(mod, what = "density", data = dms_log, breaks = 30)
abline(v=gmm_thresh, col='red', lwd=2, lty=2)
dev.off()

### get gene mat
ds_qt_all_l = c()
ds_qt_all_h = c()
for (ct in common_ct){
ct_pos_rna = which(rna_list==ct)
if (length(ct_pos_rna)>1){
rna_mat_tmp = rowMeans(ds_qt[,ct_pos_rna])
}else{rna_mat_tmp = ds_qt[,ct_pos_rna]}
ds_qt_all_l = c(ds_qt_all_l, rna_mat_tmp[rna_l])
ds_qt_all_h = c(ds_qt_all_h, rna_mat_tmp[rna_h])
}

ds_qt_all_log_l = log(ds_qt_all_l+smallnum)
ds_qt_all_log_h = log(ds_qt_all_h+smallnum)


### match high low gene means
dmslog_noMean_l = ds_qt_all_log_l #- mean(c(ds_qt_all_log_l, ds_qt_all_log_h))
dmslog_noMean_h = ds_qt_all_log_h #- mean(c(ds_qt_all_log_l, ds_qt_all_log_h))


### get state coverage
sp_all_l = c()
sp_all_h = c()
for (i in c(0:23)){
ss = read.table(paste('HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.S',i,'.bed', sep=''), header=F)
sp = (ss[,-c(1:3)])
sp_i_l = c()
sp_i_h = c()
for (ct in common_ct){
ct_pos_sp = which(sp_list==ct)
### average replicates
if (length(ct_pos_sp)>1){
sp_mat_tmp_l = rowMeans(sp[rna_l,ct_pos_sp])
sp_mat_tmp_h = rowMeans(sp[rna_h,ct_pos_sp])
}else{
sp_mat_tmp_l = sp[rna_l,ct_pos_sp]
sp_mat_tmp_h = sp[rna_h,ct_pos_sp]
}
sp_i_l = c(sp_i_l, sp_mat_tmp_l)
sp_i_h = c(sp_i_h, sp_mat_tmp_h)
}
print(i)
sp_all_l = cbind(sp_all_l, sp_i_l)
sp_all_h = cbind(sp_all_h, sp_i_h)
}

colnames(sp_all_l) = c(0:23)
colnames(sp_all_h) = c(0:23)
sp_all_log_l = sp_all_l#log(sp_all_l+smallnum_sp)
sp_all_log_h = sp_all_h#log(sp_all_h+smallnum_sp)


sp_all_dist_l = c()
sp_all_dist_h = c()
for (i in c(0:23)){
ss = read.table(paste('HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.S',i,'.bed', sep=''), header=F)
sp = (ss[,-c(1:3)])
sp_i_l = c()
sp_i_h = c()
for (ct in common_ct){
ct_pos_sp = which(sp_list==ct)
### average replicates
if (length(ct_pos_sp)>1){
sp_mat_tmp_l = rowMeans(sp[rna_l,ct_pos_sp])
sp_mat_tmp_h = rowMeans(sp[rna_h,ct_pos_sp])
}else{
sp_mat_tmp_l = sp[rna_l,ct_pos_sp]
sp_mat_tmp_h = sp[rna_h,ct_pos_sp]
}
sp_i_l = c(sp_i_l, sp_mat_tmp_l)
sp_i_h = c(sp_i_h, sp_mat_tmp_h)
}
print(i)
sp_all_dist_l = cbind(sp_all_dist_l, sp_i_l)
sp_all_dist_h = cbind(sp_all_dist_h, sp_i_h)
}

colnames(sp_all_dist_l) = c(0:23)
colnames(sp_all_dist_h) = c(0:23)

sp_all_dist_log_l = log(sp_all_dist_l+smallnum_sp)
sp_all_dist_log_h = log(sp_all_dist_h+smallnum_sp)

### prepare pca
sp_all_PD_log_l = cbind(sp_all_log_l, sp_all_dist_log_l)
sp_all_PD_log_h = cbind(sp_all_log_h, sp_all_dist_log_h)


### dimension reducation
pca_all <- prcomp(rbind(sp_all_PD_log_l,sp_all_PD_log_h)[,-c(1,25)], center = T,scale. = T)
exp_var = summary(pca_all)$importance[2,]
exp_var_sum = 0
k = 0
for (varexp_tmp in exp_var){
k = k+1
exp_var_sum = exp_var_sum+varexp_tmp
if (exp_var_sum>0.95){
pcn = k
break
}
}
### PCA lm
used_idtry = c(ds_qt_all_log_l, ds_qt_all_log_h)>(-1000000)

lm_all = lm(c(dmslog_noMean_l,dmslog_noMean_h)[used_idtry]~pca_all$x[used_idtry,c(1:pcn)])
Bpca_all = pca_all$rotation[,c(1:pcn)] %*% lm_all$coefficients[-1]
#lmvar_all = lmvar(c(dmslog_noMean_l,dmslog_noMean_h), X_mu = pca_all$x[,1:pcn], X_sigma = pca_all$x[,1:pcn], intercept_mu=F)
#Bpca_lmvar_all = pca_all$rotation[,1:pcn] %*% lmvar_all$coefficients_mu


pred = lm_all$fitted.value #+ mean(c(ds_qt_all_log_l, ds_qt_all_log_h))
#pca_all$x[,1:pcn] %*% lm_all$coefficients + mean(c(ds_qt_all_log_l, ds_qt_all_log_h))
R2(c(dmslog_noMean_l,dmslog_noMean_h)[used_idtry], pred)
0.5891845
0.558521
0.4530207
0.5144291


set.seed(2019)
dmslog_noMean = c(dmslog_noMean_l,dmslog_noMean_h)[used_idtry]
ds_qt_all_log = c(ds_qt_all_log_l, ds_qt_all_log_h)[used_idtry]
used_id = sample(length(pred), 50000)
plot_lim = plot_lim_all#c(min(c(dmslog_noMean+mean(ds_qt_all_log), pred)), max(c(dmslog_noMean+mean(ds_qt_all_log), pred)))


#+mean(ds_qt_all_log)
png('pred_obs.tpm1.wg.png')
heatscatter(as.vector(pred)[used_id], as.vector(dmslog_noMean)[used_id], xlim=plot_lim, ylim=plot_lim, main=R2(c(dmslog_noMean_l,dmslog_noMean_h)[used_idtry], pred))
abline(0,1, col='red', lwd=2)
abline(h=0, col='black', lwd=2)
abline(v=0, col='black', lwd=2)
dev.off()


### get coe
#Bpca_all = pca_all$rotation[,1:pcn] %*% lm_all$coefficients[-1] 
Bpca_all = pca_all$rotation[,c(1:pcn)] %*% lm_all$coefficients[-1]

eRP_mat_human = cbind(Bpca_all[1:23],Bpca_all[24:46])
eRP_mat_human = rbind(c(0,0), eRP_mat_human)
colnames(eRP_mat_human) = c('P','D')
rownames(eRP_mat_human) = 0:23
write.table(eRP_mat_human, 'statep_rna_coe_heatmap.human.all.wg.txt', quote=F, col.names=T, row.names=T, sep='\t')
plot_lim_P = max(abs(eRP_mat_human[,1]))
plot_lim_D = max(abs(eRP_mat_human[,2]))
plot_lim_PD = max(abs(eRP_mat_human))
pdf('statep_rna_coe_heatmap.human.all.P.pdf', width=3)
rank = c(9,2,3,1,7,4,15,14,13,5,11,22,6,18,8,24,12,10,20,16,23,17,19,21)
breaksList = seq(-plot_lim_P, plot_lim_P, by = 0.001)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList))
col_breaks = c(seq(0, 2000,length=33))
pheatmap(cbind(eRP_mat_human[rank,1],eRP_mat_human[rank,1]), color=my_colorbar, breaks = breaksList, cluster_cols = FALSE,cluster_rows=F, clustering_method = 'average',annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()
pdf('statep_rna_coe_heatmap.human.all.D.pdf', width=3)
rank = c(9,2,3,1,7,4,15,14,13,5,11,22,6,18,8,24,12,10,20,16,23,17,19,21)
breaksList = seq(-plot_lim_D, plot_lim_D, by = 0.001)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList))
col_breaks = c(seq(0, 2000,length=33))
pheatmap(cbind(eRP_mat_human[rank,2],eRP_mat_human[rank,2]), color=my_colorbar, breaks = breaksList, cluster_cols = FALSE,cluster_rows=F, clustering_method = 'average',annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()



