pvaluefunc <- function(mt)
{
  m <- t(matrix(unlist(mt), 2))
  # colnames(m) <- c("Normal", "Tumor")
  # if(length(m[m<5]) != 0) {
  pv <- fisher.test(as.matrix(m))$p.value
  #  }else {
  #pv <- chisq.test(as.matrix(m))$p.value
  # }
  return(pv)
}