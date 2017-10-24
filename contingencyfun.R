contingencyfun <- function(memb)
{
  memb[,1] <- gsub("[[:punct:]]", "-", memb[,1])
  input <- cbind(apply(memb[1], 1, function(x) ifelse(grepl("-G1", x), "Group1", "Group2")), memb) #check input is from shiny?
  
  mytable <- table(input[,4], input[,1])
  mytable <- as.data.frame.matrix(mytable)
  df <- unmatrix(mytable,byrow=T)
  return(df)
}