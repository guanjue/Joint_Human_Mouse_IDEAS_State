d1 = read.table('tpmFeb21_v3.tab', sep='\t', header=T)
d2 = read.table('annotation3.idsort.od.txt', sep='\t', header=F)

d2_id = apply(d2,1, function(x) unlist(strsplit(x[1], '[.]'))[1])

shared_d2_id = is.element(d2_id, d1[,1])

d2a = d2[shared_d2_id,]


write.table(d2a, 'annotation3.idsort.txt', quote=F, sep='\t', col.names=F, row.names=F)
