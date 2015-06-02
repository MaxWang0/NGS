qt <- read.table("/net/zmf5/cb/19/yxw124430/Review_of_Work/GeneFpk/gene_fpk_matrix_with_gi.txt", header = T, sep = "\t")                           ### read in read_depth matrix
gene.index <- qt[1:nrow(qt)-1,1]                
qt <- qt[,2:126]
qt2 <- as.vector(as.matrix(qt))                                                                 ### transform to vector from frame
dt <- c()
row <- c()                                                                                      ### define empty vector
merge <- c()                                                                                    ### define empty vector 
for (i in 1:(nrow(qt)-1)){
	row <- as.vector(as.matrix(qt[i,]))
	parameter <- lm(row ~ as.vector(as.matrix(qt[nrow(qt),]))+0)
	slope <- as.vector(as.matrix(coef(parameter)))
	r <- round(cor(row, as.vector(as.matrix(qt[nrow(qt),]))),digits = 2)
	vector <- c(slope,r)
	merge <- rbind(merge, vector)
}

merge <- data.frame(merge)
merge <- cbind(gene.index, merge)
r2 <- merge[,3]*merge[,3]
merge <- cbind(merge, r2)
for (i in 2:ncol(merge)){
	for (j in 1:nrow(merge)){
	merge[,i][j] <- round(merge[,i][j],digits=2)                                                      ###gene affinity
	}
}
write.table(merge, "nosim.uni.exon.affinity.bed", row.names = F, col.names = T, sep = "\t", quote = F)
