plotMeans <- function(data, data_order, type, dist) {
  
  lT2 = list()
  df = list()
  for(j in 1:length(type)) {
    target = as.vector(data_order[[j]])
    data2 = data[[j]][,match(target, data[[j]][1, ])]
  
    data3 = data2[-1,]
    group = factor(unlist(data3[1,]))
  
    df_mean = vector()
    df_var = vector()
    df_median = vector()
    lT = list()
  
    for(i in 1:length(unique(group))) {
      e1 = data3[,data3[1,] %in% i]
      if(!is.matrix(e1)) { 
        e2= e1[-1]
      } else if(is.matrix(e1)) { 
        e2 = e1[-1,]
      }
    
    if(dist[j] == "pearson" |dist[j] == "spearman") {
      e2.dist = as.matrix(as.dist(1-cor(e2,method=dist[j])))
    } else {
      e2.dist = dist(t(e2), method = dist[j]) 
    }
    
    lT[[i]] <- lowerTriangle(e2.dist, diag= F)
    df_mean[i] <- mean(lowerTriangle(e2.dist), na.rm = T)
    df_median[i] <- median(lowerTriangle(e2.dist), na.rm = T)
    df_var[i] <- var(lowerTriangle(e2.dist), na.rm = T)
  }
  
  df[[j]] <- data.frame(mean = df_mean, median= df_median, var = df_var)
  lT2[[j]] = lT
  }
  
  par(mfcol= c(2,length(type)))
  for(k in 1:length(type)){
    plot(df[[k]]$mean, df[[k]]$var, xlab = "Cluster Means", ylab= "Cluster Variance", main = type[k])
    text(df[[k]]$mean, df[[k]]$var, row.names(df[[k]]), cex=0.8, pos=4, col="red")
  
    boxplot(lT2[[k]], main = "Cluster distribution")
  }
  #dev.off()
  

}
