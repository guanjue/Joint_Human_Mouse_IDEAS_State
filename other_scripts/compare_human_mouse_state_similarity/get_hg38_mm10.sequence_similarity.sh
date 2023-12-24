cd /Users/guanjuexiang/Documents/projects/analysis/0813_human_mouse_state_compare_heatmap
# get fasta sequence for hg38 and mm10
bedtools getfasta -fi /Users/guanjuexiang/Documents/projects/analysis/04_05_2022_coe/coe_analysis/hg38.fa -bed hg38.gene.GATA1.matched_ct.state.bed > hg38.gene.GATA1.matched_ct.state.bed.fa
bedtools getfasta -fi /Users/guanjuexiang/Downloads/Snapshot_test/mm10.fa -bed mm10.gene.Gata1.matched_ct.state.bed > mm10.gene.Gata1.matched_ct.state.bed.fa

# open R
library(Biostrings)

# read fasta sequence
GATA1_hg38_fa = readDNAStringSet("hg38.gene.GATA1.matched_ct.state.bed.fa")
Gata1_mm10_fa = readDNAStringSet("mm10.gene.Gata1.matched_ct.state.bed.fa")

# get similarity matrix
similarity_mat = matrix(0, nrow = length(GATA1_hg38_fa), ncol = length(Gata1_mm10_fa))

# get similarity score
for (i in 1:length(GATA1_hg38_fa)){
    print(i)
  for (j in 1:length(Gata1_mm10_fa)){
    similarity_mat[i,j] = pairwiseAlignment(GATA1_hg38_fa[i], Gata1_mm10_fa[j], type = "global-local", substitutionMatrix = "BLOSUM62", gapOpening = 10, gapExtension = 1, scoreOnly = TRUE)
  }
}

# save similarity matrix
write.table(t(similarity_mat), "GATA1_hg38_Gata1_mm10_DNA_sequence_similarity.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# plot heatmap
library(pheatmap)
pdf("GATA1_hg38_Gata1_mm10_DNA_sequence_similarity_heatmap.pdf")
pheatmap(t(similarity_mat), cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE, color = colorRampPalette(c("white", "black"))(100))
dev.off()

