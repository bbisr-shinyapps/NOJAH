sil_cpt <- function(data_use, opt_k, res, dist, sil_cp, cols = cols, cpt_data, cpt_method, cpt_max) {
  
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
  rownames(sk2) = colnames(dt) # m0$order
  
  sil.order <- as.numeric(rownames(sortSilhouette(sk2)))
  
  m0.1 <- m0[match(sil.order, m0$order),]
  
  sil_width <- sk2[,3]
  sil_width_ordered <- sil_width[sil.order]
  
  m0.1 <- cbind(m0.1,sil_width_ordered )
  neg <- list()
  neg3 <- list()
  count_cluster <- list()
  
  for(i in 1:opt_k) {
    count_cluster[[i]] <- nrow(m0.1[m0.1$Cluster == i,])
    if(count_cluster[[i]] <= 5){
        stop("Number of samples are too few. Please choose another method!")
    } else { 
        neg[[i]] <- cptGenes(x= m0.1[m0.1$Cluster == i,5], cpt_data = cpt_data, cpt_method= cpt_method, cpt_max = cpt_max)
    }  
  }
  
  neg2 <- do.call(rbind, neg)
  
  m0.1 <- cbind(m0.1, neg2)
  m0.1$order2 <- 1:nrow(m0.1)
  
  
  for(i in 1:opt_k) {
    if(sil_cp[[i]] != 0) {
    neg3[[i]] <-  which( m0.1[m0.1$Cluster == i, "changepoints"] > sil_cp[[i]]) ### for each cluster find samples in CP=1
    } else if(sil_cp[[i]] == 0) {
      neg3[[i]] <- NA
    } 
  }

  neg4 <- list(neg3[[1]])
 
  for(j in 2:opt_k) {
    neg4[[2]] <- neg3[[j]] + count_cluster[[j-1]]
  }
  #neg_sil_index <- which(m0.1[,"changepoints"] > )
  neg_sil_index <- unlist(neg4)
  neg_sil_index <- neg_sil_index[!is.na(neg_sil_index)]
  
  res3 = m0.1[!(m0.1$order2 %in% neg_sil_index), ]
  m <- merge(res3, cols, by = "Sample", sort = F)
  res4 <- m[,2]
  
  data_use2 = data_use[, (colnames(data_use) %in% res3$Sample)]
  data_use2 =  data_use2[, match(m$Sample, colnames(data_use2))]
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
  
  check_sk2 <<- sk2
  check_sk2.col <<- m0$colors
  check_res3 <<- res3
  
  return(list(sk2= sk2, sk3= sk3, sk2.col = m0.1$colors, sk3.col = m$colors, core.samples = res3$Sample, k = length(unique(res3$Cluster))))
  
}
