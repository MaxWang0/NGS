qt <- read.table("gene_ODF_value.bed", header = T, sep = "\t")
odf <- as.vector(as.matrix(qt))
sort_odf <- sort(odf, decreasing = T)
len <- length(odf)
per_10 <- 0.1*len
cutoff <- sort_odf[per_10]
index <- c()
for(i in 1:125){
	index <- c(index,  sum(qt[,i]>cutoff))
}
a <- 1:125
names(a) <- index
real_index <- as.numeric(a[as.character(sort(index),decreasing = T)])
write.table(real_index, "10%_odf_index.txt", row.names = F, col.names = T, sep= "\t", quote = F)
