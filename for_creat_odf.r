qt <- read.table("/net/zmf5/cb/19/yxw124430/wy/12_12_2014/cpp/gene/new.gene.region.fpk.txt", header = T, sep = "\t")
merge <- read.table("gene_fpk_alpha.txt", header = T, sep = "\t")
ORD <- qt[,2:126]
ORD <- ORD[1:nrow(qt)-1, ]
qt <- as.vector(as.matrix(qt[nrow(qt),]))
zvalue <- c()
ERD <- c()
expected <- c()
MRD <- as.numeric(qt[2:126])
MRDt <- rep(MRD, nrow(ORD))
MRD <- matrix(MRDt, nrow= 125, ncol= nrow(ORD))
MRD <- t(MRD)
alpha <- c()
for(i in 1:125){ 
	alpha <- cbind(alpha, merge[,2])
}
expected <- alpha * MRD
zvalue <- round((ORD-expected)/sqrt(expected), digits = 2)
zvalue_ORD_expected <- cbind(zvalue, ORD, expected)
zvalue_ORD_expected <- zvalue_ORD_expected[!is.na(zvalue[,1])&!is.infinite(zvalue[,1]),]
nz <- zvalue_ORD_expected[,1:125]
nORD <- zvalue_ORD_expected[,126:250]
nexpected <- zvalue_ORD_expected[,251:375]
nz <- data.frame(nz)
nexpected <- data.frame(nexpected)
ord <- data.frame(nORD)
write.table(ord, "3_exon_ORD.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(nz, "3_exon_zvalue.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(nexpected, "3_exon_ERD.bed", row.names = F, col.names = T, sep = "\t", quote  = F)
zvalue <- read.table("3_exon_zvalue.bed", header = T, sep = "\t")
ERD <- read.table("3_exon_ERD.bed", header = T, sep = "\t")
Cg <- c()
Cs <- c()
Csg <- c()
odf <- c()
Cg <- apply(nz, 1, sd)
Cs <- apply(nz, 2, sd)

write.table(Cg, "Cg.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(Cs, "Cs.bed", row.names = F, col.names = T,sep = "\t", quote = F)
meanCg <- mean(Cg)



Cg <- as.vector(as.matrix(Cg))
Cs <- as.vector(as.matrix(Cs))
Cgt <- rep(Cg, ncol(nz))
Cst <- rep(Cs, nrow(nz))
Css <- matrix(Cst, nrow= 125, ncol= nrow(nz))
Css <- t(Css)
Cgg <- c()

for (i in 1:125) {
	Cgg <- cbind(Cgg, Cg)
}

odf <- (Cgg * Css)/meanCg

write.table(odf, "gene_ODF_value.bed", row.names = F, col.names = T, sep = "\t", quote = F)
