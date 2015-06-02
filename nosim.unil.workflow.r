qt <- read.table("/net/zmf5/cb/19/yxw124430/Review_of_Work/GeneFpk/gene_fpk_matrix_with_gi.txt", header = T, sep = "\t")                           ### read in read_depth matrix
gene.index <- read.table("/net/zmf5/cb/19/yxw124430/wy/conifer_v0.2.2/summer.topic/read_depth/rpk_RD/merged.sorted.by.gene.start.zxp_2011_10_26_0-FH_1700_16.read.depth.bed", header = F, sep = "\t")  ### read in first-column gene.index
sample.index <- read.table("/net/zmf5/cb/19/yxw124430/wy/conifer_v0.2.2/summer.topic/read_depth/unilinear_model/new2.read_depth_matrix.bed", header = F, sep = "\t")

gene.index <- qt[,1]
gene.index <- gene.index[1: length(gene.index)-1]
qt2 <- as.vector(as.matrix(qt))                                                                 ### transform to vector from frame
dt <- c()
row <- c()                                                                                      ### define empty vector
merge <- c()                                                                                    ### define empty vector 
qt <- qt[,2:ncol(qt)]
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
write.table(merge, "nosim.uni.gene.affinity.bed", row.names = F, col.names = T, sep = "\t", quote = F)
for (i in 1:(nrow(qt)-1)){
	dt <- rbind(dt,as.vector(as.matrix(qt[i,])))
}
dt <- cbind(merge[,4],data.frame(dt))                                                                 ###add r2
dt <- cbind(merge[,2],data.frame(dt))                                                                 ###add slope
dt <- cbind(merge[,1],data.frame(dt))                                                                 ###add gene
dt3 <- dt[dt[,3]>=0.7&!is.na(dt[,3]),]                                                                                ###add r2 cutoff
write.table(dt3, "nosim.uni.r2.cutoff.read_depth_matrix_with_gene_name.bed", row.names = F, col.names = T, sep = "\t", quote = F)
row <- c()
c <- c()
row1 <- c()
row0 <- c()
row2 <- c()
row3 <- c()
CN <- c()
ORD <- c()
expected <- c()
variance <- c()
cn0 <- c()
cn1 <- c()
cn2 <- c()
cn3 <- c()
cn4 <- c()
cn5 <- c()
cn6 <- c()
cn7 <- c()
cn8 <- c()
cn9 <- c()
CN1 <- c()
CN2 <- c()
CN3 <- c()
CN4 <- c()
CN5 <- c()
CN6 <- c()
CN7 <- c()
CN8 <- c()
CN9 <- c()
for (i in 1:nrow(dt3)){
	for (j in 4:ncol(dt3)){
		expected[j-3] <- qt[nrow(qt),(j-3)]*as.vector(as.matrix(dt3[i,][2]))
		c[j-3] <- round((as.vector(as.matrix(dt3[i,][j]))-expected[j-3])/sqrt(qt[nrow(qt),(j-3)]*as.vector(as.matrix(dt3[i,][2]))),digits = 2)
	}
	row <- rbind(row,c)
	row1 <- rbind(row1,expected)
}
zvalue <- data.frame(row)
 ERD <- data.frame(row1)
write.table(zvalue, "nosim.uni.z_value.bed", row.names = F, col.names =T, sep = "\t", quote = F)
write.table(data.frame(row1), "nosim.uni.expected_value.bed", row.names= F, col.names = T, sep = "\t", quote = F) 
Cg <- c()
Cs <- c()
Csg <- c()
odf <- c()
row2 <- c()
cn0 <- c()
cn1 <- c()                                                                                                                                                      
cn2 <- c()
cn3 <- c()                                                                                                                                                        
cn4 <- c()                                                                                                                                                        
cn5 <- c()                                                                                                                                                        
cn6 <- c()                                                                                                                                                        
cn7 <- c()                                                                                                                                                        
cn8 <- c()
cn9 <- c()
CN0 <- c()
CN1 <- c()
CN2 <- c()
CN3 <- c()
CN4 <- c()
CN5 <- c()
CN6 <- c()
CN7 <- c()
CN8 <- c()
CN9 <- c()
denominator <- c()
iteration <- c()
variance <- c()

for (i in 1:nrow(zvalue)){
	Cg[i] <- round(sd(as.vector(as.matrix(zvalue[i,]))), digits = 2)
}
for (j in 1:ncol(zvalue)){
	Cs[j] <- round(sd(as.vector(as.matrix(zvalue[,j]))), digits = 2) 
}
meanCg <- mean(Cg)                                                                                 
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

for (i in 1:nrow(zvalue)){
	for (j in 1:ncol(zvalue)){
		Csg[j] <- round((Cg[i] * Cs[j])/meanCg, digits = 2)
		variance[j] <- round(Csg[j] * sqrt(ERD[i,j]), digits = 2)
		P0 = (1/(variance[j]*p0))*exp(-(((dt3[,j+3][i]-ERD[,j][i]/10)/(variance[j]*sqrt[1]))^2)/2)*0.000634
		P1 = (1/(variance[j]*p1))*exp(-(((dt3[,j+3][i]-ERD[,j][i]/2)/(variance[j]*sqrt[2]))^2)/2)*0.00211
		P2 = (1/(variance[j]*p2))*exp(-(((dt3[,j+3][i]-ERD[,j][i])/(variance[j]*sqrt[3]))^2)/2)*0.996
		P3 = (1/(variance[j]*p3))*exp(-(((dt3[,j+3][i]-ERD[,j][i]*1.5)/(variance[j]*sqrt[4]))^2)/2)*0.000538
		P4 = (1/(variance[j]*p4))*exp(-(((dt3[,j+3][i]-ERD[,j][i]*2)/(variance[j]*sqrt[5]))^2)/2)*0.000668
		P5 = (1/(variance[j]*p5))*exp(-(((dt3[,j+3][i]-ERD[,j][i]*2.5)/(variance[j]*sqrt[6]))^2)/2)*0.0000357
		P6 = (1/(variance[j]*p6))*exp(-(((dt3[,j+3][i]-ERD[,j][i]*3)/(variance[j]*sqrt[7]))^2)/2)*0.00000752
		P7 = (1/(variance[j]*p7))*exp(-(((dt3[,j+3][i]-ERD[,j][i]*3.5)/(variance[j]*sqrt[8]))^2)/2)*0.00000139
		P8 = (1/(variance[j]*p8))*exp(-(((dt3[,j+3][i]-ERD[,j][i]*4)/(variance[j]*sqrt[9]))^2)/2)*0.000000361
		P9 = (1/(variance[j]*p9))*exp(-(((dt3[,j+3][i]-ERD[,j][i]*4.5)/(variance[j]*sqrt[10]))^2)/2)*0.0000000437
		sum = P0 + P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9

		cn0[j] = P0 / sum
		cn1[j] = P1 / sum
		cn2[j] = P2 / sum
		cn3[j] = P3 / sum
		cn4[j] = P4 / sum
		cn5[j] = P5 / sum
		cn6[j] = P6 / sum
		cn7[j] = P7 / sum
		cn8[j] = P8 / sum
		cn9[j] = P9 / sum
	
	}
	odf <- rbind(odf,Csg)
	row2 <- rbind(row2, variance)
	CN0 <- rbind(CN0, cn0)
	CN1 <- rbind(CN1, cn1) 
	CN2 <- rbind(CN2, cn2)
	CN3 <- rbind(CN3, cn3)
	CN4 <- rbind(CN4, cn4)
	CN5 <- rbind(CN5, cn5)
	CN6 <- rbind(CN6, cn6)
	CN7 <- rbind(CN7, cn7)
	CN8 <- rbind(CN8, cn8)
	CN9 <- rbind(CN9, cn9)
	
}
frame <- data.frame(odf)
dev <- data.frame(row2)
CN0 <- data.frame(CN0)
CN1 <- data.frame(CN1)
CN2 <- data.frame(CN2)
CN3 <- data.frame(CN3)
CN4 <- data.frame(CN4)
CN5 <- data.frame(CN5)
CN6 <- data.frame(CN6)
CN7 <- data.frame(CN7)
CN8 <- data.frame(CN8)
CN9 <- data.frame(CN9)
write.table(frame, "nosim.uni.ODF_value.bed", row.names = F, col.names = T, sep = "\t", quote= F)
write.table(dev, "nosim.uni.variance_value.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(CN0, "nosim.uni.h_value_for_CN0.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(CN1, "nosim.uni.h_value_for_CN1.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(CN2, "nosim.uni.h_value_for_CN2.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(CN3, "nosim.uni.h_value_for_CN3.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(CN4, "nosim.uni.h_value_for_CN4.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(CN5, "nosim.uni.h_value_for_CN5.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(CN6, "nosim.uni.h_value_for_CN6.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(CN7, "nosim.uni.h_value_for_CN7.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(CN8, "nosim.uni.h_value_for_CN8.bed", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(CN9, "nosim.uni.h_value_for_CN9.bed", row.names = F, col.names = T, sep = "\t", quote = F)

row <- c()
for (i in 1:ncol(CN0)){
	for (j in 1:nrow(CN0)){
		if (CN0[,i][j] >= 0.65){
			row <- rbind(row, c("CN0", j, i, CN0[,i][j]))
		}
	}
}
for (i in 1:ncol(CN1)){
	for (j in 1:nrow(CN1)){
		if (CN1[,i][j] >= 0.65){
			row <- rbind(row, c("CN1", j, i, CN1[,i][j]))
		}
	}
}
for (i in 1:ncol(CN3)){
	for (j in 1:nrow(CN3)){
		if (CN3[,i][j] >= 0.65){
			row <- rbind(row, c("CN3", j, i, CN3[,i][j]))
		}
	}
}
for (i in 1:ncol(CN4)){
	for (j in 1:nrow(CN4)){
		if (CN4[,i][j] >= 0.65){
			row <- rbind(row, c("CN4", j, i, CN4[,i][j]))
		}
	}
}
for (i in 1:ncol(CN5)){
	for (j in 1:nrow(CN5)){
		if (CN5[,i][j] >= 0.65){
			row <- rbind(row, c("CN5", j, i, CN5[,i][j]))
		}
	}
}
for (i in 1:ncol(CN6)){
	for (j in 1:nrow(CN6)){
		if (CN6[,i][j] >= 0.65){
			row <- rbind(row, c("CN6", j, i, CN6[,i][j]))
		}
	}
}
for (i in 1:ncol(CN7)){
	for (j in 1:nrow(CN7)){
		if (CN7[,i][j] >= 0.65){
			row <- rbind(row, c("CN7", j, i, CN7[,i][j]))
		}
	}
}
for (i in 1:ncol(CN8)){
	for (j in 1:nrow(CN8)){
		if (CN8[,i][j] >= 0.65){
			row <- rbind(row, c("CN8", j, i, CN8[,i][j]))
		}
	}
}
for (i in 1:ncol(CN9)){
	for (j in 1:nrow(CN9)){
		if (CN9[,i][j] >= 0.65){
			row <- rbind(row, c("CN9", j, i, CN9[,i][j]))
		}
	}
}
write.table(row,"nosim.uni.002.call_for_Bayesian.bed", sep = "\t", row.names =F, col.names=T, quote= F)

ORDmatrix <- dt3
ERDmatrix <- ERD
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
write.table(qt2, "nosim.uni.003.call_for_Bayesian.bed", row.names = F, col.names = T, sep = "\t", quote = F)
sample <- read.table("/net/zmf5/cb/19/yxw124430/wy/conifer_v0.2.2/summer.topic/read_depth/unilinear_model/new2.read_depth_matrix.bed", header = F, sep = "\t")
gene <- dt3
sam <- c()
g <- c()
qt3 <- c()
for (i in 1:nrow(qt2)){
	sam[i] <- as.vector(as.matrix(sample[1,as.numeric(qt2[i,3])]))
	g[i] <- as.vector(as.matrix(gene[qt2[i,2],1]))
}
qt3 <- cbind(qt2,g,sam)
write.table(qt3, "nosim.uni.004.call_for_Bayesian.bed", row.names = F, col.names = T, sep = "\t", quote = F)

