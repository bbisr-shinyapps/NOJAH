
coca <- function(cc, type, opt_k, coca_reps, clust, dist, coca_opt_k) {
  v.cast = list()
  for ( i in 1:length(type))
  {
  v = as.data.frame(t(cc[[i]]$output))
  v$Sample = rownames(v)
  v$value = 1
  names(v)[1] = "Cluster"
  v.cast[[i]] = cast(v,  Cluster~ Sample, value = 'value')
  v.cast[[i]][,1] = paste(type[i], 1:opt_k[i], sep = "_")
  }

  cc_datain <- do.call("rbind", v.cast)
  cc_datain[is.na(cc_datain)] <- 0
  
  cc_datain2 <- cc_datain[, -1]
  rownames(cc_datain2) = cc_datain[,1]
  
  data_in = data.matrix(cc_datain2)
  
  finalcc <- consensus_clustering(dinput=data_in, mK=10, rep=coca_reps, pI=0.8, pF= 1, cAlg="hc", dist=dist, iL=clust, fL=clust)
 
  return(list(output= finalcc[["output"]][[coca_opt_k]]$consensusClass, data= data_in, distance = finalcc[["distance"]], data_in = data_in))
}


