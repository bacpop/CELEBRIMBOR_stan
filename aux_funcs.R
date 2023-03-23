read_matrix <- function(path){
  # Read PA matrix using data.table and covert back to matrix with named rows and columns
  dt <- data.table::fread(path)
  dt <- data.table::melt(dt, id.vars = "V1")
  dt <- dt[, .(observed = sum(value), n_samples = .N), "variable"]
  data.table::setnames(dt, "variable", "gene")
  dt[, crude_freq := observed / n_samples]
  
  return(dt)
}

cel_ml <- function(y, N, c){
  dat <- list(N = N,
              y = y,
              completeness = c,
              a = 1,
              b = 1)
  fit <- mod$optimize(data = dat, init = dat$y / dat$N)
  out <- as.numeric(fit$summary("true_freq")$estimate)
  return(out)
}