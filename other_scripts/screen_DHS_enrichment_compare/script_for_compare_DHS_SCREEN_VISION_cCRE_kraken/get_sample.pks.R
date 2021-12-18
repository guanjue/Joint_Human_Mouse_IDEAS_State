args = commandArgs(trailingOnly=TRUE)


input_file = args[1]
sample_num = as.numeric(args[2])
seed = as.numeric(args[3])

d = read.table(input_file, header=F)

set.seed(seed)

ds = d[sample(dim(d)[1], sample_num),]

write.table(ds, paste(input_file, '.',seed ,'.sample.bed', sep=''), quote=F, col.names=F, row.names=F, sep='\t')


