cCRE_gene_JmetID_FS01_geneName_MID_HID_S01_TSS_GeneKM_mat_hg38_RNAcor_SFN_mat = read.table('cCRE.Gene100KB.hg38.JmetID.FS01.geneName.HID.MID.S01.TSS.GeneKMID.RNAcor.SFNID.bed', header=T, sep='\t')




#################################################
### write Gene-cCRE loop file
hg38_gene_locus = read.table('/Users/guanjuexiang/Documents/projects/analysis/0813_human_mouse_state_compare_heatmap/hg38.gene.bed', header=F)
hg38_gene_locus_TSS = hg38_gene_locus
hg38_gene_locus_TSS[hg38_gene_locus_TSS[,4]=='+',3] = hg38_gene_locus_TSS[hg38_gene_locus_TSS[,4]=='+',2]+1
hg38_gene_locus_TSS[hg38_gene_locus_TSS[,4]=='-',2] = hg38_gene_locus_TSS[hg38_gene_locus_TSS[,4]=='-',3]-1
### modify CSF1R gene locus
hg38_gene_locus_TSS[hg38_gene_locus_TSS[,5]=='CSF1R',c(2:3)] = c(150086553, 150086554)
### add TSS position
hg38_SFN_loops = cCRE_gene_JmetID_FS01_geneName_MID_HID_S01_TSS_GeneKM_mat_hg38_RNAcor_SFN_mat
hg38_SFN_loops_TSS = hg38_SFN_loops[,1:3]
k = 0
for (gene_i in unique(hg38_SFN_loops$GeneName)){
  k = k+1
  if (k%%1000==0){print(k)}
  used_n = sum(hg38_SFN_loops$GeneName==gene_i)
  TSS_gene_i = hg38_gene_locus_TSS[hg38_gene_locus_TSS[,5]==gene_i,1:3]
  hg38_SFN_loops_TSS[hg38_SFN_loops$GeneName==gene_i,] = cbind(as.data.frame(rep(TSS_gene_i[1,1], used_n)), rep(TSS_gene_i[1,2], used_n), rep(TSS_gene_i[1,3], used_n))
}
###
hg38_SFN_loops = cbind(hg38_SFN_loops[,1:5], round(hg38_SFN_loops[12]*1000), hg38_SFN_loops[,13], rep('#DEDEDE', dim(hg38_SFN_loops)[1]), hg38_SFN_loops[,c(1:3,7)], rep('.', dim(hg38_SFN_loops)[1]), hg38_SFN_loops_TSS, hg38_SFN_loops$GeneName, rep('.', dim(hg38_SFN_loops)[1]) )
### remove neg
#hg38_SFN_loops[(hg38_SFN_loops[,6]<0) | (is.na(hg38_SFN_loops[,6])),6] = 0
#hg38_SFN_loops[(is.na(hg38_SFN_loops[,6])),6] = 0
### change color

hg38_SFN_loops[hg38_SFN_loops[,7]=='SF+',8] = '#FF0000'
hg38_SFN_loops[hg38_SFN_loops[,7]=='SF',8] = '#ED585E'
hg38_SFN_loops[hg38_SFN_loops[,7]=='S+',8] = '#FF5500'
hg38_SFN_loops[hg38_SFN_loops[,7]=='S',8] = '#4FE54A'
hg38_SFN_loops[hg38_SFN_loops[,7]=='J',8] = '#FFAB00'
hg38_SFN_loops[hg38_SFN_loops[,7]=='N',8] = '#5868F2'
#hg38_SFN_loops[hg38_SFN_loops[,7]=='TSS',8] = '#000000'

hg38_SFN_loops[hg38_SFN_loops[,17]=='GATA1',]
hg38_SFN_loops[hg38_SFN_loops[,17]=='CSF1R',]
options("scipen"=100, "digits"=4)
hg38_SFN_loops$HID = apply(cbind(hg38_SFN_loops$JmetID, hg38_SFN_loops$HID), 1, function(x) paste(x[1], x[2], sep=':'))
write.table(hg38_SFN_loops, 'hg38.SFN.loop.100KB.interact', quote=F, sep='\t', col.names=F, row.names=F)

### get new loop interact file with hierachy TSS>SFJ>SF; TSS>SJ>S; TSS>J>N
### hg38
ptm = proc.time()
unique_HID = unique(hg38_SFN_loops$HID)
hg38_SFN_loops_new = hg38_SFN_loops
k = 0
for ( cCRE_i in unique_HID){
if (k%%1000==0){print(k)}
  d_cCRE_i = hg38_SFN_loops[hg38_SFN_loops$HID==cCRE_i,]
if (dim(d_cCRE_i)[1]==1){
  k = k+1
  hg38_SFN_loops_new[k,] = d_cCRE_i
} else {
  if (sum(d_cCRE_i[,7]=='TSS')>0){
    d_cCRE_i = d_cCRE_i[is.element(d_cCRE_i[,7], c('TSS')),]
  } else if (sum(d_cCRE_i[,7]=='SF+')>0) {
    d_cCRE_i = d_cCRE_i[is.element(d_cCRE_i[,7], c('SF+')),]
  } else if (sum(d_cCRE_i[,7]=='S+')>0) {
    d_cCRE_i = d_cCRE_i[is.element(d_cCRE_i[,7], c('S+')),]
  } else if (sum(d_cCRE_i[,7]=='J')>0) {
    d_cCRE_i = d_cCRE_i[is.element(d_cCRE_i[,7], c('J')),]
  }
  ###
  if (dim(d_cCRE_i)[1]==1){
    k = k+1
    hg38_SFN_loops_new[k,] = d_cCRE_i
  } else{
    k_vec = (k+1):(k+dim(d_cCRE_i)[1])
    hg38_SFN_loops_new[k_vec,] = d_cCRE_i
    k = k_vec[length(k_vec)]
  }
}
}
hg38_SFN_loops_new1 = hg38_SFN_loops_new[1:k,]
proc.time() - ptm


### write interaction passing correlation threshold
library(tidyverse)
library(grDevices)
hg38_SFN_loops_new2 = hg38_SFN_loops_new1
hg38_SFN_loops_new2[,5] = hg38_SFN_loops_new1[,6]
hg38_SFN_loops_new2[,6] = hg38_SFN_loops_new1[,5]
colnames(hg38_SFN_loops_new2)[c(6,5)] = colnames(hg38_SFN_loops_new1)[c(6,5)]
write.table(hg38_SFN_loops_new2, 'hg38.SFNJ.loop.OD.100KB.interact', quote=F, sep='\t', col.names=F, row.names=F)
### write interact file with header
bash1 = 'cat loop.header.txt > hg38.SFN.loop.interact.tmp'
bash2 = 'cat hg38.SFNJ.loop.OD.100KB.interact | awk -F \'\t\' -v OFS=\'\t\'  \'{if ($5>=0) print $0; else print $1,$2,$3,$4,0,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}\' | awk -F \'\t\' \'{if ($17=="HBA1") print $0}\' >> hg38.SFN.loop.interact.tmp && mv hg38.SFN.loop.interact.tmp hg38.SFN+.loop.HBA1.100KB.interact'
system(bash1)
system(bash2)







cCRE_gene_JmetID_FS01_geneName_MID_HID_S01_TSS_GeneKM_mat_mm10_SFN_mat = read.table('cCRE.Gene100KB.mm10.JmetID.FS01.geneName.MID.HID.S01.TSS.GeneKMID.SFNID.bed', header=T, sep='\t')



#################################################
### write Gene-cCRE loop file
mm10_gene_locus = read.table('/Users/guanjuexiang/Documents/projects/analysis/0813_human_mouse_state_compare_heatmap/mm10.gene.bed', header=F)
mm10_gene_locus_TSS = mm10_gene_locus
mm10_gene_locus_TSS[mm10_gene_locus_TSS[,4]=='+',3] = mm10_gene_locus_TSS[mm10_gene_locus_TSS[,4]=='+',2]+1
mm10_gene_locus_TSS[mm10_gene_locus_TSS[,4]=='-',2] = mm10_gene_locus_TSS[mm10_gene_locus_TSS[,4]=='-',3]-1
### add TSS position
mm10_SFN_loops = cCRE_gene_JmetID_FS01_geneName_MID_HID_S01_TSS_GeneKM_mat_mm10_SFN_mat
mm10_SFN_loops_TSS = mm10_SFN_loops[,1:3]
k = 0
for (gene_i in unique(mm10_SFN_loops$GeneName)){
  k = k+1
  if (k%%1000==0){print(k)}
  used_n = sum(mm10_SFN_loops$GeneName==gene_i)
  TSS_gene_i = mm10_gene_locus_TSS[toupper(mm10_gene_locus_TSS[,5])==gene_i,1:3]
  mm10_SFN_loops_TSS[mm10_SFN_loops$GeneName==gene_i,] = cbind(as.data.frame(rep(TSS_gene_i[1,1], used_n)), rep(TSS_gene_i[1,2], used_n), rep(TSS_gene_i[1,3], used_n))
}
###
mm10_SFN_loops = cbind(mm10_SFN_loops[,1:5], round(rep(1, dim(mm10_SFN_loops)[1])*1000), mm10_SFN_loops[,12], rep('#DEDEDE', dim(mm10_SFN_loops)[1]), mm10_SFN_loops[,c(1:3,7)], rep('.', dim(mm10_SFN_loops)[1]), mm10_SFN_loops_TSS, mm10_SFN_loops$GeneName, rep('.', dim(mm10_SFN_loops)[1]) )
### remove neg
mm10_SFN_loops[mm10_SFN_loops[,7]=='SF+',8] = '#FF0000'
mm10_SFN_loops[mm10_SFN_loops[,7]=='SF',8] = '#ED585E'
mm10_SFN_loops[mm10_SFN_loops[,7]=='S+',8] = '#FF5500'
mm10_SFN_loops[mm10_SFN_loops[,7]=='S',8] = '#4FE54A'
mm10_SFN_loops[mm10_SFN_loops[,7]=='J',8] = '#FFAB00'
mm10_SFN_loops[mm10_SFN_loops[,7]=='N',8] = '#5868F2'

mm10_SFN_loops[mm10_SFN_loops[,17]=='GATA1',]
mm10_SFN_loops[mm10_SFN_loops[,17]=='CSF1R',]
options("scipen"=100, "digits"=4)
mm10_SFN_loops$HID = apply(cbind(mm10_SFN_loops$JmetID, mm10_SFN_loops$MID), 1, function(x) paste(x[1], x[2], sep=':'))
write.table(mm10_SFN_loops, 'mm10.SFN.loop.100KB.interact', quote=F, sep='\t', col.names=F, row.names=F)

### get new loop interact file with hierachy TSS>SFJ>SF; TSS>SJ>S; TSS>J>N
### mm10
ptm = proc.time()
unique_MID = unique(mm10_SFN_loops$HID)
mm10_SFN_loops_new = mm10_SFN_loops
k = 0
for ( cCRE_i in unique_MID){
if (k%%1000==0){print(k)}
  d_cCRE_i = mm10_SFN_loops[mm10_SFN_loops$HID==cCRE_i,]
if (dim(d_cCRE_i)[1]==1){
  k = k+1
  mm10_SFN_loops_new[k,] = d_cCRE_i
} else {
  if (sum(d_cCRE_i[,7]=='TSS')>0){
    d_cCRE_i = d_cCRE_i[is.element(d_cCRE_i[,7], c('TSS')),]
  } else if (sum(d_cCRE_i[,7]=='SF+')>0) {
    d_cCRE_i = d_cCRE_i[is.element(d_cCRE_i[,7], c('SF+')),]
  } else if (sum(d_cCRE_i[,7]=='S+')>0) {
    d_cCRE_i = d_cCRE_i[is.element(d_cCRE_i[,7], c('S+')),]
  } else if (sum(d_cCRE_i[,7]=='J')>0) {
    d_cCRE_i = d_cCRE_i[is.element(d_cCRE_i[,7], c('J')),]
  }
  ###
  if (dim(d_cCRE_i)[1]==1){
    k = k+1
    mm10_SFN_loops_new[k,] = d_cCRE_i
  } else{
    k_vec = (k+1):(k+dim(d_cCRE_i)[1])
    mm10_SFN_loops_new[k_vec,] = d_cCRE_i
    k = k_vec[length(k_vec)]
  }
}
}
mm10_SFN_loops_new1 = mm10_SFN_loops_new[1:k,]
proc.time() - ptm

### write interaction passing correlation threshold
library(tidyverse)
library(grDevices)
mm10_SFN_loops_new2 = mm10_SFN_loops_new1
mm10_SFN_loops_new2[,5] = mm10_SFN_loops_new1[,6]
mm10_SFN_loops_new2[,6] = mm10_SFN_loops_new1[,5]
colnames(mm10_SFN_loops_new2)[c(6,5)] = colnames(mm10_SFN_loops_new1)[c(6,5)]
mm10_SFN_loops_new2[,12] = mm10_SFN_loops_new2[,19]
mm10_SFN_loops_new2 = mm10_SFN_loops_new2[,c(1:18)]
write.table(mm10_SFN_loops_new2, 'mm10.SFNJ.loop.OD.100KB.interact', quote=F, sep='\t', col.names=F, row.names=F)
### write interact file with header
bash1 = 'cat loop.mm10.header.txt > mm10.SFN.loop.interact.tmp'
bash2 = 'cat mm10.SFNJ.loop.OD.100KB.interact | awk -F \'\t\' -v OFS=\'\t\' \'{if ($5>=0 || $7=="TSS") print $0; else print $1,$2,$3,$4,0,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}\' | sort -u | awk -F \'\t\' \'{if ($17=="HBA1") print $0}\' >> mm10.SFN.loop.interact.tmp && mv mm10.SFN.loop.interact.tmp mm10.SFN+.loop.HBA1.100KB.interact'
system(bash1)
system(bash2)



