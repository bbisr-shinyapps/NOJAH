col_color <- function(col1, col.groups, number.col.groups, col.groups.name, col = c("darkblue", "grey", "orange", "yellow", "purple", "darkgreen", "hotpink", "brown", "darkorchid2", "maroon"))
  {
  
  if(number.col.groups==1) { 
    cell <- c(rep(col.groups.name, number.col.groups))
    cc1 <- rep(col1[50], length(cell))
    cc1[1:table(col.groups)[[1]]] <- col[1]
  } else if(number.col.groups==2) {
    cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]))
    cc1 <- rep(col1[50], length(cell))
    cc1[1:table(col.groups)[[1]]] <- col[1]
    cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- col[2]
  } else if(number.col.groups==3) {
    cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]])) 
    cc1 <- rep(col1[50], length(cell))
    cc1[1:table(col.groups)[[1]]] <- col[1]
    cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- col[2]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- col[3]
  } else if(number.col.groups==4) {
    cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]])) 
    cc1 <- rep(col1[50], length(cell))
    cc1[1:table(col.groups)[[1]]] <- col[1]
    cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- col[2]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- col[3]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- col[4]
  } else if(number.col.groups==5) {
    cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]])) 
    cc1 <- rep(col1[50], length(cell))
    cc1[1:table(col.groups)[[1]]] <- col[1]
    cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- col[2]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- col[3]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- col[4]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- col[5]
  } else if(number.col.groups==6) {
    cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]])) 
    cc1 <- rep(col1[50], length(cell))
    cc1[1:table(col.groups)[[1]]] <- col[1]
    cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- col[2]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- col[3]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- col[4]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- col[5]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- col[6]
  } else if(number.col.groups==7) {
    cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]])) 
    cc1 <- rep(col1[50], length(cell))
    cc1[1:table(col.groups)[[1]]] <- col[1]
    cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- col[2]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- col[3]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- col[4]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- col[5]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- col[6]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- col[7]
  } else if(number.col.groups==8) {
    cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]])) 
    cc1 <- rep(col1[50], length(cell))
    cc1[1:table(col.groups)[[1]]] <- col[1]
    cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- col[2]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- col[3]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- col[4]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- col[5]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- col[7]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- col[8]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- col[9]
  } else if(number.col.groups==9) {
    cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]]), rep(col.groups.name[9], table(col.groups)[[9]])) 
    cc1 <- rep(col1[50], length(cell))
    cc1[1:table(col.groups)[[1]]] <- col[1]
    cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- col[2]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- col[3]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- col[4]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- col[5]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- col[6]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- col[7]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- col[8]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+1:table(col.groups)[[9]]] <- col[9]
  } else if(number.col.groups==10) {
    cell <- c(rep(col.groups.name[1], table(col.groups)[[1]]), rep(col.groups.name[2],table(col.groups)[[2]]), rep(col.groups.name[3],table(col.groups)[[3]]), rep(col.groups.name[4],table(col.groups)[[4]]), rep(col.groups.name[5],table(col.groups)[[5]]), rep(col.groups.name[6], table(col.groups)[[6]]), rep(col.groups.name[7], table(col.groups)[[7]]), rep(col.groups.name[8], table(col.groups)[[8]]), rep(col.groups.name[9], table(col.groups)[[9]]), rep(col.groups.name[10], table(col.groups)[[10]])) 
    cc1 <- rep(col1[50], length(cell))
    cc1[1:table(col.groups)[[1]]] <- col[1]
    cc1[table(col.groups)[[1]]+1:table(col.groups)[[2]]] <- col[2]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+1:table(col.groups)[[3]]] <- col[3]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+1:table(col.groups)[[4]]] <- col[4]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+1:table(col.groups)[[5]]] <- col[5]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+1:table(col.groups)[[6]]] <- col[6]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+1:table(col.groups)[[7]]] <- col[7]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+1:table(col.groups)[[8]]] <- col[8]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+1:table(col.groups)[[9]]] <- col[9]
    cc1[table(col.groups)[[1]]+table(col.groups)[[2]]+table(col.groups)[[3]]+table(col.groups)[[4]]+table(col.groups)[[5]]+table(col.groups)[[6]]+table(col.groups)[[7]]+table(col.groups)[[8]]+table(col.groups)[[9]]+1:table(col.groups)[[10]]] <- col[10]
  }
  return(cc1)
  
}