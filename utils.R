chop <- function(x, position){
  if(position == 1){
    out <- as.numeric(substr(x = x, 
                             start = as.numeric(gregexpr(pattern = "\\[", text = x)) + 1,
                             stop = as.numeric(gregexpr(pattern = ",", text = x)) - 1))
  }else if(position == 2){
    out <- as.numeric(substr(x = x, 
                             start = as.numeric(gregexpr(pattern = ",", text = x)) + 1,
                             stop = as.numeric(gregexpr(pattern = "\\]", text = x)) - 1))
  }
  return(out)
}