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

setwd('/gpfs/group/yzz2/default/legacy/group/projects/vision_joint_human_mouse/coefficients_human_with_cCRE_with_corfilter/')

sp_list = c('B', 'B', 'CD34', 'CD34', 'CD4', 'CD4', 'CD8', 'CD8', 'CLP', 'CLP', 'CMP', 'CMP', 'EOS', 'ERY', 'ERY', 'GMP', 'GMP', 'H1', 'H1', 'H2', 'H2', 'LMPP', 'LMPP', 'LSK', 'LSK', 'MK', 'MK', 'MONc', 'MONc', 'MONp', 'MONp', 'MPP', 'MPP', 'NEU', 'NEU', 'NK', 'NK', 'PBMC', 'PBMC')
rna_list = c('LSK', 'LSK', 'ERY', 'ERY', 'CD4', 'CD4', 'CD8', 'CD8', 'B', 'B', 'CMP', 'CMP', 'MONp', 'MONp', 'NEU', 'NEU', 'MONc', 'MONc', 'GMP', 'GMP', 'CFUE', 'NK', 'NK', 'MK', 'MK', 'CLP', 'CLP', 'MPP', 'MPP', 'EOS', 'EOS')
common_ct = intersect(sp_list, rna_list)
smallnum = 1
cor_thresh = 0.5
pcv_var_used = 0.95
set.seed(2019)

### get RNA
d = read.table('../correlation/HumanVISION_RNAseq_hg38_genes_tpm.idsort.protein_coding.txt', header=F)
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
ds_qt_all = c()

for (ct in common_ct){
ct_pos_rna = which(rna_list==ct)
if (length(ct_pos_rna)>1){
rna_mat_tmp = rowMeans(ds_qt[,ct_pos_rna])
}else{rna_mat_tmp = ds_qt[,ct_pos_rna]}
ds_qt_all = c(ds_qt_all, rna_mat_tmp)
}

ds_qt_all_log = log(ds_qt_all+smallnum)

### get state coverage
sp_all = c()
for (i in c(0:23)){
ss = read.table(paste('../coefficients_human_with_cCRE/HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.Nkbupdownexp.S',i,'.bed', sep=''), header=F)
sp = (ss[,-c(1:3)])
sp_i = c()
for (ct in common_ct){
ct_pos_sp = which(sp_list==ct)
### average replicates
if (length(ct_pos_sp)>1){
sp_mat_tmp = rowMeans(sp[,ct_pos_sp])
}else{
sp_mat_tmp = sp[,ct_pos_sp]
}
sp_i = c(sp_i, sp_mat_tmp)
}
print(i)
sp_all = cbind(sp_all, sp_i)
}

colnames(sp_all) = c(0:23)
sp_all_log = log(sp_all+smallnum)


sp_all_dist = c()
for (i in c(0:23)){
ss = read.table(paste('../coefficients_human_with_cCRE/HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.S',i,'.bed', sep=''), header=F)
sp = (ss[,-c(1:3)])
sp_i = c()
for (ct in common_ct){
ct_pos_sp = which(sp_list==ct)
### average replicates
if (length(ct_pos_sp)>1){
sp_mat_tmp = rowMeans(sp[,ct_pos_sp])
}else{
sp_mat_tmp = sp[,ct_pos_sp]
}
sp_i = c(sp_i, sp_mat_tmp)
}
print(i)
sp_all_dist = cbind(sp_all_dist, sp_i)
}

colnames(sp_all_dist) = c(0:23)

sp_all_dist_log = log(sp_all_dist+smallnum)

### prepare pca
sp_all_PD_log = cbind(sp_all_log, sp_all_dist_log)

### match high low gene means
dmslog_noMean = ds_qt_all_log - mean(ds_qt_all_log)


### dimension reducation
pca_all <- prcomp(sp_all_PD_log[,-c(1,25)], center = F,scale. = F)
exp_var = summary(pca_all)$importance[2,]
exp_var_sum = 0
k = 0
for (varexp_tmp in exp_var){
k = k+1
exp_var_sum = exp_var_sum+varexp_tmp
if (exp_var_sum>pcv_var_used){
pcn = k
break
}
}
### PCA lm
lm_all = lm(dmslog_noMean~pca_all$x[,1:pcn]-1)
Bpca_all = pca_all$rotation[,1:pcn] %*% lm_all$coefficients
lmvar_all = lmvar(dmslog_noMean, X_mu = pca_all$x[,1:pcn], X_sigma = pca_all$x[,1:pcn], intercept_mu=F)
Bpca_lmvar_all = pca_all$rotation[,1:pcn] %*% lmvar_all$coefficients_mu

pred = pca_all$x[,1:pcn] %*% lm_all$coefficients + mean(ds_qt_all_log)
R2(ds_qt_all_log, pred)
0.4077193

eRP_mat_human = cbind(Bpca_lmvar_all[1:23],Bpca_lmvar_all[24:46])
eRP_mat_human = rbind(c(0,0), eRP_mat_human)
colnames(eRP_mat_human) = c('P','D')
rownames(eRP_mat_human) = 0:23
write.table(eRP_mat_human, 'statep_rna_coe_heatmap.human.all.txt', quote=F, col.names=T, row.names=T, sep='\t')

used_id = sample(length(pred), 10000)
plot_lim = c(min(c(dmslog_noMean+mean(ds_qt_all_log), pred)), max(c(dmslog_noMean+mean(ds_qt_all_log), pred)))
png('pred_obs.tpm1.png')
heatscatter(as.vector(dmslog_noMean+mean(ds_qt_all_log))[used_id], as.vector(pred)[used_id], xlim=plot_lim, ylim=plot_lim)
dev.off()


### filter by cor
### get state coe predict mat
ss00 = read.table(paste('bp24_v2_cCRE_max1.cCRE.M.withid.S0.mat.txt', sep=''), header=F)

sp_sig_pred_ct = matrix(0, nrow=dim(ss00)[1], ncol=dim(ss00)[2])
for (i in c(0:23)){
print(i)
ss = read.table(paste('bp24_v2_cCRE_max1.cCRE.M.withid.S',i,'.mat.txt', sep=''), header=F)
sp = (ss[,-c(1:4)])
sp_sig = log(sp+smallnum)*eRP_mat_human[i+1,2]
sp_sig_pred_ct = sp_sig_pred_ct + sp_sig
}

sp_sig_pred_ct_common_ct = c()
for (ct in common_ct){
ct_pos_sp = which(sp_list==ct)
### average replicates
if (length(ct_pos_sp)>1){
sp_mat_tmp = rowMeans(sp_sig_pred_ct[,ct_pos_sp])
}else{
sp_mat_tmp = sp_sig_pred_ct[,ct_pos_sp]
}
sp_sig_pred_ct_common_ct = cbind(sp_sig_pred_ct_common_ct, sp_mat_tmp)
}

### get gene mat
ds_qt_all_mat = c()

for (ct in common_ct){
ct_pos_rna = which(rna_list==ct)
if (length(ct_pos_rna)>1){
rna_mat_tmp = rowMeans(ds_qt[,ct_pos_rna])
}else{rna_mat_tmp = ds_qt[,ct_pos_rna]}
ds_qt_all_mat = cbind(ds_qt_all_mat, rna_mat_tmp)
}

ds_qt_all_mat_log = log(ds_qt_all_mat+smallnum)


### read gene with used id
gene_withccre_id = read.table('HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.withccreid.bed', header=F)[,5]
gene_withccre_id_str = as.character(gene_withccre_id)

cor_vec_all_select = matrix('.', nrow=length(gene_withccre_id), ncol=1)

for (i in 1:length(gene_withccre_id)){
#for (i in 1:6){
if (i%%1000==0){
	print(i)
}
ds_qt_all_mat_log_i = ds_qt_all_mat_log[i,]
state_id = as.numeric(unlist(strsplit(gene_withccre_id_str[i], ',')))
sp_sig_pred_ct_common_ct_i = sp_sig_pred_ct_common_ct[state_id,]
if (length(state_id)>1){
cor_vec_i = apply(sp_sig_pred_ct_common_ct_i, 1, function(x) cor(x, ds_qt_all_mat_log_i))
} else if (length(state_id)==1){
cor_vec_i = cor(sp_sig_pred_ct_common_ct_i, ds_qt_all_mat_log_i)
}
cor_vec_i_select = state_id[(cor_vec_i>=cor_thresh) & (!is.na(cor_vec_i))]
cor_vec_i_select_str = paste(cor_vec_i_select, collapse=',')
cor_vec_all_select[i,1] = cor_vec_i_select_str
}

### redo sp coverage after cor filtering
ss = read.table(paste('bp24_v2_cCRE_max1.cCRE.M.withid.S0.mat.txt', sep=''), header=F)
ss_all = (ss[,-c(1:4)])
for (i in c(1:23)){
print(i)
ss = read.table(paste('bp24_v2_cCRE_max1.cCRE.M.withid.S',i,'.mat.txt', sep=''), header=F)
sp = (ss[,-c(1:4)])
ss_all = cbind(ss_all, sp)
}

sp_corfiltered = matrix(0, nrow=dim(cor_vec_all_select)[1], ncol=dim(ss_all)[2])

for (i in 1:dim(cor_vec_all_select)[1]){
if (i%%100==0){
	print(i)
}
used_row_i = as.numeric(unlist(strsplit(cor_vec_all_select[i,1], ',')))
if (length(used_row_i)>0){
sp_corfiltered[i,] = colSums(ss_all[used_row_i,])
}else{
sp_corfiltered[i,] = rep(0, dim(sp_corfiltered)[2])
}
}


sp_corfiltered_all = matrix(0, nrow=dim(sp_corfiltered)[1]*length(common_ct), ncol=24)
k=0
for (ct in common_ct){
print(ct)
ct_pos_sp = which(sp_list==ct)
sp_mat_tmp = matrix(0, nrow=dim(sp_corfiltered)[1], ncol=24 )
for (i in ct_pos_sp){
sp_mat_tmp = sp_mat_tmp+sp_corfiltered[,seq(1,39*24, by=39)+i-1]
}
sp_mat_tmp = sp_mat_tmp/length(ct_pos_sp)
sp_corfiltered_all[seq(1,dim(sp_corfiltered)[1])+(k)*dim(sp_corfiltered)[1],] = sp_mat_tmp
k = k+1
}

colnames(sp_corfiltered_all) = c(0:23)
sp_corfiltered_all_D_log = log(sp_corfiltered_all+smallnum)

sp_all_PD_corfilter_log = cbind(sp_all_log, sp_corfiltered_all_D_log)


### dimension reducation rm 0 state info -c(1,25)
pca_all <- prcomp(sp_all_PD_corfilter_log[,-c(1,25)], center = F,scale. = F)
exp_var = summary(pca_all)$importance[2,]
exp_var_sum = 0
k = 0
for (varexp_tmp in exp_var){
k = k+1
exp_var_sum = exp_var_sum+varexp_tmp
if (exp_var_sum>pcv_var_used){
pcn = k
break
}
}
### PCA lm
lm_all = lm(dmslog_noMean~pca_all$x[,1:pcn]-1)
Bpca_all = pca_all$rotation[,1:pcn] %*% lm_all$coefficients
#lmvar_all = lmvar(dmslog_noMean, X_mu = pca_all$x[,1:pcn], X_sigma = pca_all$x[,1:pcn], intercept_mu=F)
#Bpca_lmvar_all = pca_all$rotation[,1:pcn] %*% lmvar_all$coefficients_mu

pred = pca_all$x[,1:pcn] %*% lm_all$coefficients + mean(ds_qt_all_log)
R2(ds_qt_all_log, pred)
0.4116009

used_id = sample(length(pred), 10000)
plot_lim = c(min(c(dmslog_noMean+mean(ds_qt_all_log), pred)), max(c(dmslog_noMean+mean(ds_qt_all_log), pred)))
png('pred_obs.tpm2.png')
heatscatter(as.vector(dmslog_noMean+mean(ds_qt_all_log))[used_id], as.vector(pred)[used_id], xlim=plot_lim, ylim=plot_lim)
dev.off()



eRP_mat_human = cbind(Bpca_lmvar_all[1:23],Bpca_lmvar_all[24:46])
eRP_mat_human = rbind(c(0,0), eRP_mat_human)
colnames(eRP_mat_human) = c('P','D')
rownames(eRP_mat_human) = 0:23
write.table(eRP_mat_human, 'statep_rna_coe_heatmap.human.all.ccre.withcorfilter.txt', quote=F, col.names=T, row.names=T, sep='\t')


#############################
#############################
#############################Round2
### filter by cor
### get state coe predict mat
sp_sig_pred_ct = matrix(0, nrow=dim(ss00)[1], ncol=dim(ss00)[2])
for (i in c(0:23)){
print(i)
ss = read.table(paste('bp24_v2_cCRE_max1.cCRE.M.withid.S',i,'.mat.txt', sep=''), header=F)
sp = (ss[,-c(1:4)])
sp_sig = log(sp+smallnum)*eRP_mat_human[i+1,2]
sp_sig_pred_ct = sp_sig_pred_ct + sp_sig
}

sp_sig_pred_ct_common_ct = c()
for (ct in common_ct){
ct_pos_sp = which(sp_list==ct)
### average replicates
if (length(ct_pos_sp)>1){
sp_mat_tmp = rowMeans(sp_sig_pred_ct[,ct_pos_sp])
}else{
sp_mat_tmp = sp_sig_pred_ct[,ct_pos_sp]
}
sp_sig_pred_ct_common_ct = cbind(sp_sig_pred_ct_common_ct, sp_mat_tmp)
}

### get gene mat
ds_qt_all_mat = c()

for (ct in common_ct){
ct_pos_rna = which(rna_list==ct)
if (length(ct_pos_rna)>1){
rna_mat_tmp = rowMeans(ds_qt[,ct_pos_rna])
}else{rna_mat_tmp = ds_qt[,ct_pos_rna]}
ds_qt_all_mat = cbind(ds_qt_all_mat, rna_mat_tmp)
}

ds_qt_all_mat_log = log(ds_qt_all_mat+smallnum)


### read gene with used id
gene_withccre_id = read.table('HumanVISION_RNAseq_hg38_gene.idsort.protein_coding.NHkbupdownexp.withccreid.bed', header=F)[,5]
gene_withccre_id_str = as.character(gene_withccre_id)

cor_vec_all_select = matrix('.', nrow=length(gene_withccre_id), ncol=1)

for (i in 1:length(gene_withccre_id)){
#for (i in 1:6){
if (i%%1000==0){
	print(i)
}
ds_qt_all_mat_log_i = ds_qt_all_mat_log[i,]
state_id = as.numeric(unlist(strsplit(gene_withccre_id_str[i], ',')))
sp_sig_pred_ct_common_ct_i = sp_sig_pred_ct_common_ct[state_id,]
if (length(state_id)>1){
cor_vec_i = apply(sp_sig_pred_ct_common_ct_i, 1, function(x) cor(x, ds_qt_all_mat_log_i))
} else if (length(state_id)==1){
cor_vec_i = cor(sp_sig_pred_ct_common_ct_i, ds_qt_all_mat_log_i)
}
cor_vec_i_select = state_id[(cor_vec_i>=cor_thresh) & (!is.na(cor_vec_i))]
cor_vec_i_select_str = paste(cor_vec_i_select, collapse=',')
cor_vec_all_select[i,1] = cor_vec_i_select_str
}

### redo sp coverage after cor filtering
ss = read.table(paste('bp24_v2_cCRE_max1.cCRE.M.withid.S0.mat.txt', sep=''), header=F)
ss_all = (ss[,-c(1:4)])
for (i in c(1:23)){
print(i)
ss = read.table(paste('bp24_v2_cCRE_max1.cCRE.M.withid.S',i,'.mat.txt', sep=''), header=F)
sp = (ss[,-c(1:4)])
ss_all = cbind(ss_all, sp)
}

sp_corfiltered = matrix(0, nrow=dim(cor_vec_all_select)[1], ncol=dim(ss_all)[2])

for (i in 1:dim(cor_vec_all_select)[1]){
if (i%%100==0){
	print(i)
}
used_row_i = as.numeric(unlist(strsplit(cor_vec_all_select[i,1], ',')))
if (length(used_row_i)>0){
sp_corfiltered[i,] = colSums(ss_all[used_row_i,])
}else{
sp_corfiltered[i,] = rep(0, dim(sp_corfiltered)[2])
}
}


sp_corfiltered_all = matrix(0, nrow=dim(sp_corfiltered)[1]*length(common_ct), ncol=24)
k=0
for (ct in common_ct){
print(ct)
ct_pos_sp = which(sp_list==ct)
sp_mat_tmp = matrix(0, nrow=dim(sp_corfiltered)[1], ncol=24 )
for (i in ct_pos_sp){
sp_mat_tmp = sp_mat_tmp+sp_corfiltered[,seq(1,39*24, by=39)+i-1]
}
sp_mat_tmp = sp_mat_tmp/length(ct_pos_sp)
sp_corfiltered_all[seq(1,dim(sp_corfiltered)[1])+(k)*dim(sp_corfiltered)[1],] = sp_mat_tmp
k = k+1
}

colnames(sp_corfiltered_all) = c(0:23)
sp_corfiltered_all_D_log = log(sp_corfiltered_all+smallnum)

sp_all_PD_corfilter_log = cbind(sp_all_log, sp_corfiltered_all_D_log)


### dimension reducation rm 0 state info -c(1,25)
pca_all <- prcomp(sp_all_PD_corfilter_log[,-c(1,25)], center = F,scale. = F)
exp_var = summary(pca_all)$importance[2,]
exp_var_sum = 0
k = 0
for (varexp_tmp in exp_var){
k = k+1
exp_var_sum = exp_var_sum+varexp_tmp
if (exp_var_sum>pcv_var_used){
pcn = k
break
}
}
### PCA lm
lm_all = lm(dmslog_noMean~pca_all$x[,1:pcn]-1)
Bpca_all = pca_all$rotation[,1:pcn] %*% lm_all$coefficients
lmvar_all = lmvar(dmslog_noMean, X_mu = pca_all$x[,1:pcn], X_sigma = pca_all$x[,1:pcn], intercept_mu=F)
Bpca_lmvar_all = pca_all$rotation[,1:pcn] %*% lmvar_all$coefficients_mu

pred = pca_all$x[,1:pcn] %*% lm_all$coefficients + mean(ds_qt_all_log)
R2(ds_qt_all_log, pred)
0.4191341


eRP_mat_human = cbind(Bpca_lmvar_all[1:23],Bpca_lmvar_all[24:46])
eRP_mat_human = rbind(c(0,0), eRP_mat_human)
colnames(eRP_mat_human) = c('P','D')
rownames(eRP_mat_human) = 0:23
write.table(eRP_mat_human, 'statep_rna_coe_heatmap.human.all.ccre.withcorfilter.txt', quote=F, col.names=T, row.names=T, sep='\t')

plot_lim_P = max(abs(eRP_mat_human[,1]))
plot_lim_D = max(abs(eRP_mat_human[,2]))
plot_lim_PD = max(abs(eRP_mat_human))

pdf('statep_rna_coe_heatmap.human.all.ccre.withcorfilter.pdf', width=3)
rank = c(9,2,3,1,7,4,15,14,13,5,11,22,6,18,8,24,12,10,20,16,23,17,19,21)
breaksList = seq(-plot_lim_PD, plot_lim_PD, by = 0.001)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList))
col_breaks = c(seq(0, 2000,length=33))
pheatmap(eRP_mat_human[rank,], color=my_colorbar, breaks = breaksList, cluster_cols = FALSE,cluster_rows=F, clustering_method = 'average',annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()



pdf('statep_rna_coe_heatmap.human.all.P.ccre.withcorfilter.pdf', width=3)
rank = c(9,2,3,1,7,4,15,14,13,5,11,22,6,18,8,24,12,10,20,16,23,17,19,21)
breaksList = seq(-plot_lim_P, plot_lim_P, by = 0.001)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList))
col_breaks = c(seq(0, 2000,length=33))
pheatmap(cbind(eRP_mat_human[rank,1],eRP_mat_human[rank,1]), color=my_colorbar, breaks = breaksList, cluster_cols = FALSE,cluster_rows=F, clustering_method = 'average',annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()

pdf('statep_rna_coe_heatmap.human.all.D.ccre.withcorfilter.pdf', width=3)
rank = c(9,2,3,1,7,4,15,14,13,5,11,22,6,18,8,24,12,10,20,16,23,17,19,21)
breaksList = seq(-plot_lim_D, plot_lim_D, by = 0.001)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList))
col_breaks = c(seq(0, 2000,length=33))
pheatmap(cbind(eRP_mat_human[rank,2],eRP_mat_human[rank,2]), color=my_colorbar, breaks = breaksList, cluster_cols = FALSE,cluster_rows=F, clustering_method = 'average',annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()

write.table(sp_all_PD_corfilter_log, 'sp_all_PD_corfilter_log.human.txt', quote=F, col.names=F, row.names=F, sep='\t')
write.table(dmslog_noMean, 'dmslog_noMean.human.txt', quote=F, col.names=F, row.names=F, sep='\t')



