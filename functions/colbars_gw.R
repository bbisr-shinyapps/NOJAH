colbars_gw <- function(df2, colscheme = 2){
  names.of.cols = colnames(df2)
  names.of.cols = names.of.cols[-1]
  n= nrow(df2)
  df3 <- NULL
  colours <- list()
  uni <- list()
  
  for(i in 2:ncol(df2))
  {  df3 <- c(df3,as.character(df2[,i]))
     uni[[i]] <- unique(df2[,i])
  }
  
  
  if(any(duplicated(as.character(unlist(uni))))) { # FALSE
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
  
  } else {
    
    no.of.levels2 <- list()
    level.names2 <- list()
    colours[[1]] <- c("darkblue", "grey", "orange", "yellow", "purple" , "darkgreen",  "hotpink", "brown", "darkorchid2", "maroon")
    colours[[2]] <- c("red", "green", "dodgerblue2", "violetred4", "yellow2", "plum4", "palevioletred4")
    colours[[3]] <- c("hotpink", "orange", "deepskyblue", "firebrickred", "brown4", "orchid4", "pink2") #topo.colors(10)
    colours[[4]] <- rainbow(10)
    colours[[5]] <- cm.colors(10)
    
    colors_used <- list()
    
   df3 <- data.matrix(df2[,-1])
   for(i in 1:ncol(df3)) {
    
    no.of.levels2[[i]] <- length(unique(df3[,i]))
    level.names2[[i]] <- sort(unique(df3[,i]))
    
     
    
    for(j in 1:length(df3[,i]))
    { 
      for(k in 1:no.of.levels2[[i]])
      {
        if(df3[j, i] == level.names2[[i]][k])
        { df3[j, i] = colours[[i]][k]
       
        break;
        }
        else {
          df3[j, i]=  df3[j, i] }
        
      }
      colors_used[[i]] <- colours[[i]][1:k]
      
    }
   }
   
   
   df4 <- matrix(df3, nrow = n)
   colnames(df4) <- names.of.cols
   return(list(data= df4, nlevels= unlist(no.of.levels2), levelname = unlist(level.names2), colors = unlist(colors_used)))
  } 
    
  
  
}