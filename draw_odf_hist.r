exon_count_ODF <- read.table("/net/zmf5/cb/19/yxw124430/wy/8_8_2014/exon/exon_ODF_value.bed", header = T, sep = "\t")
exon_fpk_ODF <- read.table("/net/zmf5/cb/19/yxw124430/wy/10_8/simulation/new_ev_best_cutoff/fpk_exon_odf.bed", header = T, sep = "\t")
gene_count_ODF <- read.table("/net/zmf5/cb/19/yxw124430/wy/8_8_2014/Daimgene0.3/nor2.ODF_value.bed", header = T, sep = "\t")
gene_fpk_ODF <- read.table("/net/zmf5/cb/19/yxw124430/wy/10_8/simulation/fpk_gene_odf.bed", header = T, sep = "\t")
gene_fpk <- read.table("/net/zmf5/cb/19/yxw124430/wy/10_8/gene.region.fpk.txt", header = T, sep = "\t")
exon_fpk <- read.table("/net/zmf5/cb/19/yxw124430/wy/10_8/exome125.region.fpk.txt", header = T, sep = "\t")




exon_count_ODF <- as.vector(as.matrix(exon_count_ODF))
exon_fpk_ODF <- as.vector(as.matrix(exon_fpk_ODF))
gene_count_ODF <- as.vector(as.matrix(gene_count_ODF))
gene_fpk_ODF <- as.vector(as.matrix(gene_fpk_ODF))

pdf("odf_histogram.pdf")

hist(gene_count_ODF, xlim = c(0,10), xlab = "gene_count_odf", main = "gene_count_odf", breaks = 1000)

hist(gene_fpk_ODF, xlim = c(0,1), xlab = "gene_fpk_odf", main = "gene_fpk_odf", breaks = 100)

hist(exon_count_ODF, xlim = c(0,10), xlab = "exon_count_odf", main = "exon_count_odf", breaks = 1000)

hist(exon_fpk_ODF, xlim = c(0,1), xlab = "exon_fpk_odf", main = "exon_fpk_odf", breaks = 2000)

dev.off()

pdf("gene_exon_length_dis.pdf")


hist(gene_fpk[,2], xlab = "gene_length", main = "Gene_length", breaks = 100)

hist(exon_fpk[,2], xlab = "exon_length", main = "Exon_length", breaks = 100)


gene_ave_len <- mean(gene_fpk[,2])

exon_ave_len <- mean(exon_fpk[,2])

ave_exon_count_ODF <- mean(exon_count_ODF)

ave_exon_fpk_ODF <- mean(exon_fpk_ODF)

ave_gene_count_ODF <- mean(gene_count_ODF)

ave_gene_fpk_ODF <- mean(gene_fpk_ODF)

total <- c(gene_ave_len, exon_ave_len, ave_exon_count_ODF, ave_exon_fpk_ODF, ave_gene_count_ODF, ave_gene_fpk_ODF)

rbind(total)
write.table(total, "avg_odf_len.txt", row.names = F, col.names = T, sep = "\t", quote = F)
