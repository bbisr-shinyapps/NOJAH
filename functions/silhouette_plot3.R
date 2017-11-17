silhouette_plot3 <- function(data_use, opt_k, res, dist, upto_width, cols = cols){
  rownames(data_use) <- NULL
  if(dist == "pearson" | dist == "spearman")
  { dt = as.dist(1-cor(data_use,method=dist))  }
  else 
  { dt = dist(t(data_use), method = dist) }
  
  #get order of samples
  res2 = as.data.frame(res)
  #res2$Sample = row.names(res2)
  res2$order = 1:nrow(res2)
  #names(res2)[1] <- "Cluster"
  
  #add colors in same order
  #m0 <- merge(res2, cols, by = "Sample", sort = F)
  res1 <- as.integer(res2$Cluster)
  
  sk2 <- silhouette(res1, dt ) # res - has the order of samples is same as in this object
  rownames(sk2) = colnames(dt)
  
  neg_sil_index <- which(sk2[, "sil_width"] < upto_width)
  
  res3 = res2[!(res2$order %in% neg_sil_index), ]
  #m <- merge(res3, cols, by = "Sample", sort = F)
  res4 <- res3$Cluster                           
  
  data_use2 = data_use[, colnames(data_use) %in% res3$Sample]
  check_data_use2 <<- data_use2
  
  if(dist == "pearson" | dist == "spearman")
  { dt4 = as.dist(1-cor(data_use2,method=dist))  }
  else 
  { dt4 = dist(t(data_use2), method = dist) }
  
  sk3   <- silhouette(res4, dt4)
  rownames(sk3) = rownames(dt4)
  
  #par(mfrow = c(1, 2))
  #plot(sk2,  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = m0$colors)
  #plot(sk3,  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = m$colors)
  
  return(list(sk2= sk2, sk3= sk3, core.samples = res3$Sample, k = length(unique(res3$Cluster))))
  
}
