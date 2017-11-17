bxplot <- function(data, data_order, type)
{
  target = as.vector(t(data_order$x))
  data2 = data[,match(target, data[1, ])]
  
  data3 = data2[-1,]
  group = factor(unlist(data3[1,]))
  
  pdf(paste(type, "_clusters.pdf", sep = ""), width = 12, height = 8)
  boxplot(data3[-1,])
  dev.off()
  
  e3 <- list()
  
  for(i in 1:length(unique(group))) {
    e1 = data3[,data3[1,] %in% i]
    if(!is.data.frame(e1)) { e2= e1[-1]}
    else { e2 = e1[-1,]}
    pdf(paste(type, "_cluster", i, ".pdf", sep= ""), width = 12, height = 8)
    boxplot(e2)
    dev.off()
    e3[[i]] <- as.numeric(unlist(e2))
  }
  
  pdf(paste(type, "_clusters_combined.pdf", sep = ""), width = 12, height = 8)
  boxplot(e3)
  dev.off()
  
}