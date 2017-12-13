consensus_clustering <- function(dinput, mK, rep, pI, pF, cAlg, dist, iL, fL) {
  
  data2 <- dinput[complete.cases(dinput),]
  data3 <- data.matrix(data2)
  
  ## GETTING GENE MEDIAN CENTER D USING PEARSON CORRLN DISTANCE FOR ConsensusClusteringPlus
  data4 <-  sweep(data3, 1, apply(data3, 1, median, na.rm = T))
  
  results = ConsensusClusterPlus(d=data4,maxK=mK, reps=rep,pItem=pI,pFeature=pF,clusterAlg=cAlg,distance=dist,innerLinkage= iL,  finalLinkage= fL, seed=1262118388.71279)
  
 return(list(output= results, data= data4, distance=dist))
  
}



