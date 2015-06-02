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
for (i in 1:nrow(nz)){
	Cg[i] <- round(sd(as.vector(as.matrix(nz[i,]))), digits = 2)
}

for (j in 1:ncol(nz)){
	Cs[j] <- round(sd(as.vector(as.matrix(nz[,j]))), digits = 2)
}

write.table(Cg, "Cg.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(Cs, "Cs.bed", row.names = F, col.names = T,sep = "\t", quote = F)
meanCg <- mean(Cg)



Cg <- as.vector(as.matrix(Cg))
Cs <- as.vector(as.matrix(Cs))
Csg <- c()
for (i in 1:nrow(nz)){
	for (j in 1:ncol(nz)){
		Csg[j] <- round((Cg[i] * Cs[j])/meanCg, digits = 2) 
	}
	odf <- rbind(odf, Csg)
	}
	frame <- data.frame(odf)
write.table(frame, "exon_ODF_value.bed", row.names = F, col.names = T, sep = "\t", quote = F)
variance <- round(frame * sqrt(ERD), digits = 2)
p0 <- round(sqrt(0.1)*sqrt(2*pi), digits = 2)
p1 <- round(sqrt(pi), digits = 2)
p2 <- round(sqrt(2*pi), digits = 2)
p3 <- round(sqrt(3*pi), digits = 2)
p4 <- round(sqrt(4*pi), digits = 2)
p5 <- round(sqrt(5*pi), digits =2)
p6 <- round(sqrt(6*pi), digits = 2)
p7 <- round(sqrt(7*pi), digits = 2)
p8 <- round(sqrt(8*pi), digits = 2)
p9 <- round(sqrt(9*pi), digits = 2)
sqrt <- c(sqrt(0.1), sqrt(0.5), 1, sqrt(1.5), sqrt(2), sqrt(2.5), sqrt(3), sqrt(3.5), 2, sqrt(4.5))
P0 = (1/(variance*p0))*exp(-(((ord-ERD/10)/(variance*sqrt[1]))^2)/2)*0.000634
P1 = (1/(variance*p1))*exp(-(((ord-ERD/2)/(variance*sqrt[2]))^2)/2)*0.00211
P2 = (1/(variance*p2))*exp(-(((ord-ERD)/(variance*sqrt[3]))^2)/2)*0.996
P3 = (1/(variance*p3))*exp(-(((ord-ERD*1.5)/(variance*sqrt[4]))^2)/2)*0.000538
P4 = (1/(variance*p4))*exp(-(((ord-ERD*2)/(variance*sqrt[5]))^2)/2)*0.000668
P5 = (1/(variance*p5))*exp(-(((ord-ERD*2.5)/(variance*sqrt[6]))^2)/2)*0.0000357
P6 = (1/(variance*p6))*exp(-(((ord-ERD*3)/(variance*sqrt[7]))^2)/2)*0.00000752
P7 = (1/(variance*p7))*exp(-(((ord-ERD*3.5)/(variance*sqrt[8]))^2)/2)*0.00000139
P8 = (1/(variance*p8))*exp(-(((ord-ERD*4)/(variance*sqrt[9]))^2)/2)*0.000000361
P9 = (1/(variance*p9))*exp(-(((ord-ERD*4.5)/(variance*sqrt[10]))^2)/2)*0.0000000437
sum = P0 + P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9

cn0 = P0 / sum
cn1 = P1 / sum
cn2 = P2 / sum
cn3 = P3 / sum
cn4 = P4 / sum
cn5 = P5 / sum
cn6 = P6 / sum
cn7 = P7 / sum
cn8 = P8 / sum
cn9 = P9 / sum

write.table(cn0, "uni.100g.1s.h_value_for_CN0.bed", row.names = F, col.names= T, sep = "\t", quote = F)
write.table(cn1, "uni.100g.1s.h_value_for_CN1.bed", row.names = F, col.names= T, sep = "\t", quote = F) 
write.table(cn2, "uni.100g.1s.h_value_for_CN2.bed", row.names = F, col.names= T, sep = "\t", quote = F) 
write.table(cn3, "uni.100g.1s.h_value_for_CN3.bed", row.names = F, col.names= T, sep = "\t", quote = F) 
write.table(cn4, "uni.100g.1s.h_value_for_CN4.bed", row.names = F, col.names= T, sep = "\t", quote = F) 
write.table(cn5, "uni.100g.1s.h_value_for_CN5.bed", row.names = F, col.names= T, sep = "\t", quote = F) 
write.table(cn6, "uni.100g.1s.h_value_for_CN6.bed", row.names = F, col.names= T, sep = "\t", quote = F) 
write.table(cn7, "uni.100g.1s.h_value_for_CN7.bed", row.names = F, col.names= T, sep = "\t", quote = F) 
write.table(cn8, "uni.100g.1s.h_value_for_CN8.bed", row.names = F, col.names= T, sep = "\t", quote = F)
write.table(cn9, "uni.100g.1s.h_value_for_CN9.bed", row.names = F, col.names= T, sep = "\t", quote = F) 

row <- c()
for (i in 1:ncol(cn0)){
	for (j in 1:nrow(cn0)){
		if (cn0[,i][j] >= 0.65){
			row <- rbind(row, c("CN0", j, i, cn0[,i][j]))
		}
	}
}
for (i in 1:ncol(cn1)){
	for (j in 1:nrow(cn1)){
		if (cn1[,i][j] >= 0.65){
			row <- rbind(row, c("CN1", j, i, cn1[,i][j]))
		}
	}
}
for (i in 1:ncol(cn3)){
	for (j in 1:nrow(cn3)){
		if (cn3[,i][j] >= 0.65){
			row <- rbind(row, c("CN3", j, i, cn3[,i][j]))
		}
	}
}
for (i in 1:ncol(cn4)){
	for (j in 1:nrow(cn4)){
		if (cn4[,i][j] >= 0.65){
			row <- rbind(row, c("CN4", j, i, cn4[,i][j]))
		}
	}
}
for (i in 1:ncol(cn5)){
	for (j in 1:nrow(cn5)){
		if (cn5[,i][j] >= 0.65){
			row <- rbind(row, c("CN5", j, i, cn5[,i][j]))
		}
	}
}
for (i in 1:ncol(cn6)){
	for (j in 1:nrow(cn6)){
		if (cn6[,i][j] >= 0.65){
			row <- rbind(row, c("CN6", j, i, cn6[,i][j]))
		}
	}
}
for (i in 1:ncol(cn7)){
	for (j in 1:nrow(cn7)){
		if (cn7[,i][j] >= 0.65){
			row <- rbind(row, c("CN7", j, i, cn7[,i][j]))
		}
	}
}
for (i in 1:ncol(cn8)){
	for (j in 1:nrow(cn8)){
		if (cn8[,i][j] >= 0.65){
			row <- rbind(row, c("cn8", j, i, cn8[,i][j]))
		}
	}
}
for (i in 1:ncol(cn9)){
	for (j in 1:nrow(cn9)){
		if (cn9[,i][j] >= 0.65){
			row <- rbind(row, c("cn9", j, i, cn9[,i][j]))
		}
	}
}
write.table(row,"uni.100g.1s.002.call_for_Bayesian.bed", sep = "\t", row.names =F, col.names=T, quote= F)

ORDmatrix <- ord
ERDmatrix <- nexpected
r2 <- c()
ORD <- c()
ERD <- c()
qt2 <- c()
for (i in 1:nrow(row)){
	ORD[i] <- as.vector(as.matrix(ORDmatrix[row[i,2],as.numeric(row[i,3])+3]))
	ERD[i] <- as.vector(as.matrix(ERDmatrix[as.numeric(row[i,2]),as.numeric(row[i,3])]))
	r2[i] <- as.vector(as.matrix(ORDmatrix[row[i,2],3]))
}
qt2 <- cbind(row,ORD,ERD,r2)
write.table(qt2, "uni.100g.1s.003.call_for_Bayesian.bed", row.names = F, col.names = T, sep = "\t", quote = F)
sample <- read.table("sample_title.bed", header = F, sep = "\t")
exon <- ord 
sam <- c()
g <- c()
qt3 <- c()
for (i in 1:nrow(qt2)){
	sam[i] <- as.vector(as.matrix(sample[1,as.numeric(qt2[i,3])]))
	g[i] <- as.vector(as.matrix(gene[qt2[i,2],1]))
}
qt3 <- cbind(qt2,g,sam)
write.table(qt3, "uni.100g.1s.004.call_for_Bayesian.bed", row.names = F, col.names = T, sep = "\t", quote = F)
B1 <- read.table("tril.uni.1000g.3s.numi.bed", header = T, sep = "\t")
n = 0
n2 = 0
n3 = 0
n4 = 0
n5 = 0
for (i in 1:nrow(B1)){ 
	for(j in 1:nrow(qt3)){
		if ( i == qt3[j,1] ){
			for (k in 1:(ncol(B1)-2)){
				if (B1[i,k] == qt3[j,3]){
					n = n + 1;
				}
			}
		}
	}
}

spec <- n/nrow(qt3)
sen <- n/((ncol(B1)-1)*nrow(B1))

