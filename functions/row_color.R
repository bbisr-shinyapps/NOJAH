row_color <- function(col1, row.groups, number.row.groups, row.groups.name, col = c("black", "gray", "hotpink", "brown1", "cyan", "maroon"))
{
  if(number.row.groups==1) { 
    cell2 <- c(rep(row.groups.name, number.row.groups))
    cc2 <- rep(col1[50], length(cell2))
    cc2[1:table(row.groups)[[1]]] <- col[1]
  } else if(number.row.groups==2) {
    cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]))
    cc2 <- rep(col1[50], length(cell2))
    cc2[1:table(row.groups)[[1]]] <- col[1]
    cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- col[2]
  } else if(number.row.groups==3) {
    cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]])) 
    cc2 <- rep(col1[50], length(cell2))
    cc2[1:table(row.groups)[[1]]] <- col[1]
    cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- col[2]
    cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- col[3]
  } else if(number.row.groups==4) {
    cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]])) 
    cc2 <- rep(col1[50], length(cell2))
    cc2[1:table(row.groups)[[1]]] <- col[1]
    cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- col[2]
    cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- col[3]
    cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- col[4]
  } else if(number.row.groups==5) {
    cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]]), rep(row.groups.name[5],table(row.groups)[[5]])) 
    cc2 <- rep(col1[50], length(cell2))
    cc2[1:table(row.groups)[[1]]] <- col[1]
    cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- col[2]
    cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- col[3]
    cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- col[4]
    cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+1:table(row.groups)[[5]]] <- col[5]
  } else if(number.row.groups==6) {
    cell2 <- c(rep(row.groups.name[1], table(row.groups)[[1]]), rep(row.groups.name[2],table(row.groups)[[2]]), rep(row.groups.name[3],table(row.groups)[[3]]), rep(row.groups.name[4],table(row.groups)[[4]]), rep(row.groups.name[5],table(row.groups)[[5]]), rep(row.groups.name[6],table(row.groups)[[6]])) 
    cc2 <- rep(col1[50], length(cell2))
    cc2[1:table(row.groups)[[1]]] <- col[1]
    cc2[table(row.groups)[[1]]+1:table(row.groups)[[2]]] <- col[2]
    cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+1:table(row.groups)[[3]]] <- col[3]
    cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+1:table(row.groups)[[4]]] <- col[4]
    cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+1:table(row.groups)[[5]]] <- col[5]
    cc2[table(row.groups)[[1]]+table(row.groups)[[2]]+table(row.groups)[[3]]+table(row.groups)[[4]]+table(row.groups)[[5]]+1:table(row.groups)[[6]]] <- col[6]
  }
  return(cc2)
}