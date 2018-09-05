silhouette_plot2 <- function(data_use, opt_k, res, dist, upto_width, cols = cols){
  rownames(data_use) <- NULL
  if(dist == "pearson" | dist == "spearman")
  { dt = as.dist(1-cor(data_use,method=dist))  }
  else 
  { dt = dist(t(data_use), method = dist) }
  
  #get order of samples
  res2 = as.data.frame(res)
  res2$Sample = row.names(res2)
  res2$order = 1:nrow(res2)
  names(res2)[1] <- "Cluster"
  
  #add colors in same order
  m0 <- merge(res2, cols, by = "Sample", sort = F)
  res1 <- m0[,2]
  
  sk2   <- silhouette(res1, dt ) # res - has the order of samples is same as in this object
  rownames(sk2) = names(dt) # m0$order
  
  sil.order <- as.numeric(rownames(sortSilhouette(sk2)))
  m0.1 <- m0[match(sil.order, m0$order),]
  
  sil_width <- sk2[,3]
  sil_width_ordered <- sil_width[sil.order]
  
  m0.1 <- cbind(m0.1,sil_width_ordered )
  neg <- list()
  
  
  for(i in 1:opt_k) {
    clust <- m0.1[m0.1$Cluster == i,]
    neg[[i]] <- ifelse(clust$sil_width < upto_width[[i]], 1, 2)
  }
  
  neg2 <- unlist(neg)
  
  m0.1 <- cbind(m0.1, neg2)
  m0.1 <- cbind(m0.1, order2 = 1:nrow(m0.1))
  
  neg_sil_index <- which(m0.1[, "neg2"] == 1)
 
  res3 = m0.1[!(m0.1$order2 %in% neg_sil_index), ]
  m <- merge(res3, cols, by = "Sample", sort = F)
  res4 <- m[,2]
  
  data_use2 <- data_use[, (colnames(data_use) %in% m$Sample)]
  data_use2 <- data_use2[, match(m$Sample, colnames(data_use2))]

  if(dist == "pearson" | dist == "spearman")
    { dt4 = as.dist(1-cor(data_use2,method=dist))  }
    else 
    { dt4 = dist(t(data_use2), method = dist) }
 
  sk3   <- silhouette(res4, dt4)
  rownames(sk3) = rownames(dt4)
  
  #par(mfrow = c(1, 2))
  #plot(sk2,  main = "Silhouette Plot of 'ALL' Samples", cex.names=0.6, max.strlen= 8, col = m0$colors)
  #plot(sk3,  main = "Silhouette Plot of 'CORE' Samples", cex.names=0.6, max.strlen= 8, col = m$colors)
  
  check_sk2 <<- sk2
  check_sk2.col <<- m0$colors
  check_res3 <<- res3
  
  return(list(sk2= sk2, sk3= sk3, sk2.col = m0.1$colors, sk3.col = m$colors.x, core.samples = res3$Sample, k = length(unique(res3$Cluster))))
  
}
