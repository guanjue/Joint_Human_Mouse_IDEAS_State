d = read.table('rnaTPMall_withcoordinates.0.txt', header=T)

ds = 2^(d[,-c(1:6)])
ds = ds-1e-3
ds[ds<0] = 0

dnew = cbind(d[,1:6],ds)
write.table(dnew, 'rnaTPM_withcoordinates.txt', quote=F, col.names=T, row.names=F, sep='\t')

