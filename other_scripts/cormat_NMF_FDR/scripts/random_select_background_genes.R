# get parameters
args = commandArgs(trailingOnly=TRUE)
working_dir = args[1]
gene_pool_file = args[2]
random_seed = as.integer(args[3])
tss_up = as.integer(args[4])
tss_down = as.integer(args[5])
n_gene = as.integer(args[6])
output_folder = args[7]

# set working directory
setwd(working_dir)

# set random seed
set.seed(random_seed)

# read in mm10_gene
mm10_gene = read.table(gene_pool_file, header = F)

# randomly pick n genes
mm10_gene_sample = mm10_gene[sample(1:nrow(mm10_gene), n_gene),]
colnames(mm10_gene_sample) = NULL

# get TSS
mm10_gene_sample_tss_pos = mm10_gene_sample[mm10_gene_sample[,4]=='+',c(1,2,2,4,5)]
mm10_gene_sample_tss_neg = mm10_gene_sample[mm10_gene_sample[,4]=='-',c(1,3,3,4,5)]
# pool together
mm10_gene_sample_tss = mm10_gene_sample
mm10_gene_sample_tss[1:nrow(mm10_gene_sample_tss_pos),] = mm10_gene_sample_tss_pos
mm10_gene_sample_tss[(nrow(mm10_gene_sample_tss_pos)+1):nrow(mm10_gene_sample_tss),] = mm10_gene_sample_tss_neg

# expand TSS
mm10_gene_sample_tss_coord = as.matrix(mm10_gene_sample_tss[,2:3])
mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='+',1] = mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='+',1] - tss_up
mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='+',2] = mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='+',2] + tss_down
mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='-',1] = mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='-',1] - tss_down
mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='-',2] = mm10_gene_sample_tss_coord[mm10_gene_sample_tss[,4]=='-',2] + tss_up
mm10_gene_sample_tss_coord[mm10_gene_sample_tss_coord<0] = 0

# get TSS
mm10_gene_sample_tss[,2:3] = mm10_gene_sample_tss_coord

# create output folder if not exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# write out TSS
for (i in 1:nrow(mm10_gene_sample_tss)) {
  write.table(mm10_gene_sample_tss[i,], paste0(output_folder, '/', 'bg.gene.', mm10_gene_sample_tss[i,5], '.bed'), col.names = F, row.names = F, quote = F, sep = '\t')
}
