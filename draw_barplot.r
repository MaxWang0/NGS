qt11 <- read.table("woitest11.summary.bed", header = F, sep = "\t")                                                                                                                                                    
qt10 <- read.table("woitest10.summary.bed", header  = F, sep = "\t")                                                                                                                                                   
qt9 <- read.table("woitest9.summary.bed", header = F, sep = "\t")                                                                                                                                                      
qt8 <- read.table("woitest8.summary.bed", header = F, sep = "\t")                                                                                                                                                      
qt5 <- read.table("woitest5.summary.bed", header = F, sep = "\t")                                                                                                                                                      qt2 <- rbind(qt9[,2], qt10[,3], qt9[,3],qt11[,3], qt8[,3],qt5[,3])                                                                                                                                                     
qt3 <- rbind(qt10[,4], qt9[,4], qt11[,4], qt8[,4], qt5[,4])                                                                                                                                                            barplot(qt2, col=c("green", "blue", "purple", "yellow", "orange", "red"), ylim = c(0,1), main = "sensitivity rate", beside = T, legend = c("uni","30%-100%tri", "20%-100%tri", "10%-100%tri", "20%-80%tri", "25%-75%tri"), args.legend = list(x = "topright", cex = .7), names.arg = c("1s","2s","3s","4s","5s","6s","7s","8s","9s","10s"))                                                                                                   
barplot(qt3, col=c("blue", "purple", "yellow", "orange", "red"), ylim = c(0,1), main = "sensitivity rate", beside = T, legend = c("uni + 30%-100%tri", "uni + 20%-100%tri", "uni + 10%-100%tri", "uni + 20%-80%tri", "uni + 25%-75%tri"), args.legend = list(x = "topright", cex = .7), names.arg = c("1s","2s","3s","4s","5s","6s","7s","8s","9s","10s")) 
