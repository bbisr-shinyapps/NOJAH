silhouette_plot <- function(data_use, opt_k, res, dist){
  rownames(data_use) <- NULL
  if(dist == "pearson" | dist == "spearman")
  { dt = as.dist(1-cor(data_use,method=dist))  }
  else 
  { dt = dist(t(data_use), method = dist) }
  
  sk2   <- silhouette(res, dt )
  plot(sk2,  main = "Silhouette Plot of ALL Samples", cex.names=0.6, max.strlen= 8)
}
