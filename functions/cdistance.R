#setwd("C:/Documents/NOJAH Images/")
#data <- read.csv("server_data.csv", header = T, sep = ",", stringsAsFactors = F, row.names=1)

cdistance <- function(data) {
data2 <- data.matrix(data)
rowidx <- seq_len(nrow(data2))

canb.dist <- function(x, j) sum((abs(x-j))/(abs(x)+abs(j)))

triangular <- combn(rowidx, 2, function(x) c(x[1], x[2], canb.dist(data2[x[1],], data2[x[2],])))

mat <- as.matrix(sparseMatrix(
  i=c(rowidx, triangular[1,], triangular[2,]), 
  j=c(rowidx, triangular[2,], triangular[1,]),
  x=c(rep(0, length(rowidx)), triangular[3,], triangular[3,])
))

colnames(mat) <- rownames(data2)
rownames(mat) <- rownames(data2)

proxy::dist(mat, method = canb.dist)

return(as.dist(mat))

}