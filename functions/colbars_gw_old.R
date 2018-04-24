colbars_gw <- function(df2, colscheme = 2){
  names.of.cols = colnames(df2)
  names.of.cols = names.of.cols[-1]
  n= nrow(df2)
  df3 <- NULL
  
  for(i in 2:ncol(df2))
  {  df3 <- c(df3,as.character(df2[,i])) }
  
  no.of.levels <- length(unique(df3))
  level.names <- sort(unique(df3))
  
  if(colscheme != 2)
  { colors <- topo.colors(no.of.levels)
  } else {
    colors <- c("darkblue", "grey", "orange", "yellow", "purple" , "darkgreen",  "hotpink", "brown", "darkorchid2", "maroon")[1:no.of.levels]
  }
  
  
  for(j in 1:length(df3))
  { 
    for(k in 1:no.of.levels)
    {
      if(df3[j] == level.names[k])
      { df3[j] =colors[k]
      break;
      }
      else {
        df3[j]= df3[j] }
    }
  }
  
  df4 <- matrix(df3, nrow = n)
  colnames(df4) <- names.of.cols
  return(list(data= df4, nlevels= no.of.levels, levelname = level.names, colors = colors))
  
}