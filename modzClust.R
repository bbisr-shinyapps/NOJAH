library(matrixStats)
library(lsr)

meanabsdev <- function(x, by )
{
  if(by == "row") meanAD <- apply(x , 1, aad)
  if(by == "col") meanAD <- apply(x , 2, aad)
  return(meanAD)
}

                 
modzClust <- function(xx, scale, zlim)
{
  if(scale != "none")
  {
    if (scale=="row") {
      if(any(length(rowMads(xx))) == 0) { zz <- (xx - rowMedians(xx)) / (rowMads(xx))} 
      else { zz <- (xx - rowMedians(xx)) / 1.253314*(meanabsdev(xx, by = "row")) }
    } else if (scale=="col") {
      if(any(length(colMads(xx))) == 0) {zz <- (xx - colMedians(xx)) / (colMads(xx)) }
      else {zz <- (xx - colMedians(xx)) / 1.253314*(meanabsdev(xx, by = "col")) }
    } else if (scale=="both") {    
        zz <- (xx - rowMedians(xx)) / rowMads(xx) #row 
        zz <- (zz - colMedians(zz)) / colMads(zz) # column scaling
   }
    zz <- pmin(pmax(zz,zlim[1]), zlim[2])
    return(list(data=zz))
  }
  #else {
  # return(list(data= as.matrix(data)))
  #}
}