library(cmdstanr)
library(ggplot2)
library(data.table)
source("utils.R")

# Number of genomes
N <- 200
# Number of genes
K <- 25
# Number of clusters
C <- 3

# Generate random completeness for each sample
completeness <- runif(N, 0.8, 1)

# Generate true frequencies at cluster level
# beta distribution parameterised to give
# nearly 0 or nearly 1 
true_freq <- matrix(nrow = C, ncol = K)
for(i in 1:C){
  true_freq[i, ] <- rbeta(K, 0.15, 0.15)
}

# Generate random cluster relative sizes
clust_rel_size <- gtools::rdirichlet(1, rep(2, 3))

# Assign samples to clusters based on sizes
clust <- sample.int(n = C, size = N, 
                    replace = TRUE, 
                    prob = clust_rel_size)

# Generates observed data from completeness and true frequencies
dat <- matrix(nrow = N, ncol = K)
for(i in 1:N){
  for(j in 1:K){
    dat[i, j] <- rbinom(1, size = 1, prob = completeness[i] * true_freq[clust[i], j])
  }
}

# List of data for stan model
data <- list(N = N,
             K = K,
             C = C,
             dat = dat,
             comp = completeness,
             clust = clust)

# Compile and run model
mod <- cmdstan_model("cluster_model.stan")
fit <- mod$sample(data = data, parallel_chains = 4)

# Get posterior samples and format them
res <- fit$draws("true_freq")
dt <- as.data.table(res)
dt[, clust := chop(variable, 1)]
dt[, gene := chop(variable, 2)]
dtall <- dt

# 95% credible intervals for estimated frequency
dt <- dt[, .(lq = quantile(value, 0.025), 
             uq = quantile(value, 0.975)), c("gene", "clust")]

# Create dataframe for the actual values
real_data <- data.frame(freq = as.vector(t(true_freq)),
                        gene = rep(1:K, C),
                        clust = rep(1:C, rep(K, C)))

# Plot estimates against actual
ggplot(data = dt, aes(x = gene, ymin = lq, ymax = uq)) +
  geom_errorbar() +
  geom_point(data = real_data, inherit.aes = FALSE,
             aes(x = gene, y = freq), col = "dodgerblue") +
  labs(y = "True frequency", x = "Gene ID") +
  facet_wrap(~clust)


# Aggregates frequencies to the population level (over all clusters)
real_all <- data.frame(freq = apply(X = true_freq, MARGIN = 2, FUN = function(x){
  return(sum(x * clust_rel_size))
}),
gene = 1:K)

# Aggregate model output to population level
dtall <- dtall[, .(gene_freq = sum(value * clust_rel_size)), c("iteration", "gene", "chain")
][, .(lq = quantile(gene_freq, 0.025), 
      uq = quantile(gene_freq, 0.975)), c("gene")]

# Plot population level actual vs estimates
ggplot(dtall, aes(x = gene, ymin = lq, ymax = uq)) + 
  geom_errorbar() +
  geom_point(data = real_all, inherit.aes = FALSE,
             aes(x = gene, y = freq), col = "dodgerblue") +
  labs(y = "True frequency", x = "Gene ID")

