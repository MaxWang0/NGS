library(Daim)
qt <- read.table("/net/zmf5/cb/19/yxw124430/Review_of_Work/GeneFpk/gene_fpk_matrix_with_gi.txt", header = T, sep = "\t")
ORD <- read.table("/net/zmf5/cb/19/yxw124430/wy/03_11_2015/3_exon_ORD_gi.bed", header = T, sep = "\t")
ERD <- read.table("/net/zmf5/cb/19/yxw124430/wy/03_11_2015/3_exon_ERD_gi.bed", header = T, sep = "\t")
ODF <- read.table("/net/zmf5/cb/19/yxw124430/wy/03_11_2015/gene_ODF_value_gi.bed", header = T, sep = "\t")
MRD <- as.vector(as.matrix(qt[nrow(qt), 2:126]))
CNcall <- function(x,y,z){  ### x ORD y ERD z ODF
	var <- round(sqrt(y), digit = 2)
	### p = (1:9/10) * pi
	### prior <- c(0.000634,0.00211,0.996,0.000538,0.000668,0.0000357,0.00000752,0.00000139,0.000000361,0.0000000437)
	### 1/z*var*p
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
	P0 = (1/(z*var*p0))*exp(-(((x-y/10)/(z*var*sqrt(0.1)))^2)/2)*0.000634
	P1 = (1/(z*var*p1))*exp(-(((x-y/2)/(z*var*sqrt(0.5)))^2)/2)*0.001
	P2 = (1/(z*var*p2))*exp(-(((x-y)/(z*var*1))^2)/2)*0.996
	P3 = (1/(z*var*p3))*exp(-(((x-y*1.5)/(z*var*sqrt(1.5)))^2)/2)*0.000538
	P4 = (1/(z*var*p4))*exp(-(((x-y*2)/(z*var*sqrt(2)))^2)/2)*0.000668
	P5 = (1/(z*var*p5))*exp(-(((x-y*2.5)/(z*var*sqrt(2.5)))^2)/2)*0.0000357
	P6 = (1/(z*var*p6))*exp(-(((x-y*3)/(z*var*sqrt(3)))^2)/2)*0.00000752
	P7 = (1/(z*var*p7))*exp(-(((x-y*3.5)/(z*var*sqrt(3.5)))^2)/2)*0.00000139
	P8 = (1/(z*var*p8))*exp(-(((x-y*4)/(z*var*sqrt(4)))^2)/2)*0.000000361
	P9 = (1/(z*var*p9))*exp(-(((x-y*4.5)/(z*var*sqrt(4.5)))^2)/2)*0.0000000437
	sum = P0 + P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9
	cn0 = round(P0 / sum, 5)
	cn1 = round(P1 / sum, 5)
	cn2 = round(P2 / sum, 5)
	cn3 = round(P3 / sum, 5)
	cn4 = round(P4 / sum, 5)
	cn5 = round(P5 / sum, 5)
	cn6 = round(P6 / sum, 5)
	cn7 = round(P7 / sum, 5)
	cn8 = round(P8 / sum, 5)
	cn9 = round(P9 / sum, 5)
	return (cn1)
	}

CNcall_3 <- function(x,y,z){  ### x ORD y ERD z ODF                                                                                                     
    var <- round(sqrt(y), digit = 2)

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
    P0 = (1/(z*var*p0))*exp(-(((x-y/10)/(z*var*sqrt(0.1)))^2)/2)*0.000634                                                                             
    P1 = (1/(z*var*p1))*exp(-(((x-y/2)/(z*var*sqrt(0.5)))^2)/2)*0.001                                                                                 
    P2 = (1/(z*var*p2))*exp(-(((x-y)/(z*var*1))^2)/2)*0.996                                                                                           
    P3 = (1/(z*var*p3))*exp(-(((x-y*1.5)/(z*var*sqrt(1.5)))^2)/2)*0.000538                                                                            
    P4 = (1/(z*var*p4))*exp(-(((x-y*2)/(z*var*sqrt(2)))^2)/2)*0.000668                                                                                
    P5 = (1/(z*var*p5))*exp(-(((x-y*2.5)/(z*var*sqrt(2.5)))^2)/2)*0.0000357                                                                               
	P6 = (1/(z*var*p6))*exp(-(((x-y*3)/(z*var*sqrt(3)))^2)/2)*0.00000752                                                                              
    P7 = (1/(z*var*p7))*exp(-(((x-y*3.5)/(z*var*sqrt(3.5)))^2)/2)*0.00000139                                                                          
    P8 = (1/(z*var*p8))*exp(-(((x-y*4)/(z*var*sqrt(4)))^2)/2)*0.000000361                                                                             
    P9 = (1/(z*var*p9))*exp(-(((x-y*4.5)/(z*var*sqrt(4.5)))^2)/2)*0.0000000437                                                                        
    sum = P0 + P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9                                                                                             
    cn0 = round(P0 / sum, 5)                                                                                                                          
    cn1 = round(P1 / sum, 5)                                                                                                                          
    cn2 = round(P2 / sum, 5)                                                                                                                          
    cn3 = round(P3 / sum, 5)                                                                                                                          
    cn4 = round(P4 / sum, 5)                                                                                                                          
    cn5 = round(P5 / sum, 5)                                                                                                                          
    cn6 = round(P6 / sum, 5)                                                                                                                          
    cn7 = round(P7 / sum, 5)                                                                                                                          
    cn8 = round(P8 / sum, 5)                                                                                                                          
    cn9 = round(P9 / sum, 5)
    
	return(cn3)
}
CNcall_4 <- function(x,y,z){  ### x ORD y ERD z ODF                                                                                                     
    var <- round(sqrt(y), digit = 2)                                                                                                                        
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
    P0 = (1/(z*var*p0))*exp(-(((x-y/10)/(z*var*sqrt(0.1)))^2)/2)*0.000634                                                                             
    P1 = (1/(z*var*p1))*exp(-(((x-y/2)/(z*var*sqrt(0.5)))^2)/2)*0.001                                                                                 
    P2 = (1/(z*var*p2))*exp(-(((x-y)/(z*var*1))^2)/2)*0.996                                                                                           
    P3 = (1/(z*var*p3))*exp(-(((x-y*1.5)/(z*var*sqrt(1.5)))^2)/2)*0.000538                                                                            
    P4 = (1/(z*var*p4))*exp(-(((x-y*2)/(z*var*sqrt(2)))^2)/2)*0.000668                                                                                
    P5 = (1/(z*var*p5))*exp(-(((x-y*2.5)/(z*var*sqrt(2.5)))^2)/2)*0.0000357                                                                               
	P6 = (1/(z*var*p6))*exp(-(((x-y*3)/(z*var*sqrt(3)))^2)/2)*0.00000752                                                                              
    P7 = (1/(z*var*p7))*exp(-(((x-y*3.5)/(z*var*sqrt(3.5)))^2)/2)*0.00000139                                                                          
    P8 = (1/(z*var*p8))*exp(-(((x-y*4)/(z*var*sqrt(4)))^2)/2)*0.000000361                                                                             
    P9 = (1/(z*var*p9))*exp(-(((x-y*4.5)/(z*var*sqrt(4.5)))^2)/2)*0.0000000437                                                                        
    sum = P0 + P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9                                                                                             
    cn0 = round(P0 / sum, 5)                                                                                                                          
    cn1 = round(P1 / sum, 5)                                                                                                                          
    cn2 = round(P2 / sum, 5)                                                                                                                          
    cn3 = round(P3 / sum, 5)                                                                                                                          
    cn4 = round(P4 / sum, 5)                                                                                                                          
    cn5 = round(P5 / sum, 5)                                                                                                                          
    cn6 = round(P6 / sum, 5)                                                                                                                          
    cn7 = round(P7 / sum, 5)                                                                                                                          
    cn8 = round(P8 / sum, 5)                                                                                                                          
    cn9 = round(P9 / sum, 5)
    
	return(cn4)
}
h1 <- CNcall(ORD[, 2:126], ERD[, 2:126], ODF[, 2:126])
h3 <- CNcall_3(ORD[, 2:126], ERD[, 2:126], ODF[, 2:126])
h4 <- CNcall_4(ORD[, 2:126], ERD[, 2:126], ODF[, 2:126])
call1 <- c()
for (i in 1: nrow(h1)){
	for(j in 1: ncol(h1)){
		if ( !is.nan(h1[i,j])){
			if (h1[i,j] >= 0.65){
			call1 <- rbind(call1, c(as.character(ORD[i,1]), "CN1", i, j, h1[i,j])) ### display gene sample and CN number
		}
	}
}
}

for (i in 1: nrow(h1)){
	for(j in 1: ncol(h1)){
		if ( !is.nan(h3[i,j])){
			if (h3[i,j] >= 0.65){ 
				call1 <- rbind(call1, c(as.character(ORD[i,1]), "CN3", i, j, h3[i,j]))
			}
		}
	}
}

for (i in 1: nrow(h1)){
	for(j in 1: ncol(h1)){ 
		if ( !is.nan(h4[i,j])){
			if (h4[i,j] >= 0.65){
				call1 <- rbind(call1, c(as.character(ORD[i,1]), "CN4", i, j, h4[i,j]))
			}
		}
	}
}
write.table(call1, "sm_detected20.bed", row.names = F, col.names = T, sep = "\t", quote = F)

