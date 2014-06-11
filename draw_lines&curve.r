qt0 <- read.table("sm_specificity16.bed", header = T, sep = "\t")                                                                                                                                                      
qt1 <- read.table("sm_specificity11.bed", header = T, sep = "\t")                                                                                                                                                      
qt2 <- read.table("sm_specificity12.bed", header = T, sep = "\t")                                                                                                                                                      
qt3 <- read.table("sm_specificity13.bed", header = T, sep = "\t")                                                                                                                                                      
qt4 <- read.table("sm_specificity14.bed", header = T, sep = "\t")                                                                                                                                                      
qt5 <- read.table("sm_specificity15.bed", header = T, sep = "\t")                                                                                                                                                      
qt8 <- read.table("sm_specificity20.bed", header = T, sep = "\t")                                                                                                                                                      
qt10 <- read.table("sm_specificity17.bed", header = T, sep = "\t")                                                                                                                                                     
qt15 <- read.table("sm_specificity19.bed", header = T, sep = "\t")                                                                                                                                                     
qt20 <- read.table("sm_specificity18.bed", header = T, sep = "\t")                                                                                                                                                     
spec <- as.numeric(c(qt0,qt1,qt2,qt3,qt4,qt5,qt8,qt10,qt15,qt20))                                                                                                                                                      
sd <- c(0,1,2,3,4,5,8,10,15,20)                                                                                                                                                                                        
pdf("spec_vs_sd.pdf")                                                                                                                                                                                                  
                                                                                                                                                                                                                       
plot(sd,spec,xlab = "C_alpha", ylab = "specificity", pch = 20)                                                                                                                                                         
xx <- seq(0,25, length = 50)                                                                                                                                                                                           
fit2 <- lm(spec~poly(sd,2,raw=TRUE))                                                                                                                                                                                   
lines(xx, predict(fit2, data.frame(sd=xx)), col="blue") 
