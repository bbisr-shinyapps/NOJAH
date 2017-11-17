zClust <- function(x, scale, zlim)
{
  if(scale != "none")
  {
    if (scale=="row") z <- t(scale(t(x)))
    if (scale=="col") z <- scale(x)
    if (scale=="both") {
      z <- t(scale(t(x))) #row scaling
      z <- scale(z) # column scaling
    }
    z <- pmin(pmax(z,zlim[1]), zlim[2])
    return(list(data=z))
  }
  #else {
  #return(list(data=as.numeric(as.matrix(data))))
  #}
}

