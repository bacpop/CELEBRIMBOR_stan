read_matrix <- function(path){
  x <- readLines(path)
  cols <- stri_split(str = x[1], regex = ",", omit_empty = TRUE)[[1]] 
  n <- length(x) - 1
  n2 <- length(cols)
  mat <- matrix(nrow = n, ncol = n2)
  rows <- c()
  
  print(paste0("0/",n))
  for(i in 1:n){
    
    if(i %% 100 == 0){
      print(paste0(i, "/", n))
    }
    temp <- stri_split(str = x[i + 1], regex = ",", omit_empty = TRUE)[[1]] 
    rows[i] <- temp[1]
    mat[i, ] <- as.numeric(temp[-1])
  }
  
  dt <- data.table(gene = cols, x = colSums(mat))
  dt[, crude_est := x / n]
  dt[, n := n]
  
  return(dt)
}