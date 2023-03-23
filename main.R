library(cmdstanr)
library(ggplot2)
library(data.table)
library(parallel)
library(magrittr)
source("aux_funcs.R")

# Read example matrix with 13466 genes and 20 genome samples
dt <- read_matrix("presence_absence_matrix.txt")

# Compile stan model
mod <- cmdstanr::cmdstan_model("pbsingle.stan")

# Generate random completeness for genome samples
set.seed(23032023)
completeness <- runif(20, 0.8, 0.99)

# Select genes with >50% crude frequency
sub_dt <- dt[crude_freq > 0.5]

# Parallel apply to get Celebrimbor maximum likelihood estimate for each gene in subset of data
sub_dt$cel_freq <- unlist(parallel::mclapply(X = sub_dt$observed, 
                                            FUN = cel_ml, 
                                            mc.cores = 8, 
                                            N = 20, 
                                            c = completeness))


# Merge all data back together
dt <- data.table::merge.data.table(dt, sub_dt, c("gene", "observed", "n_samples", "crude_freq"), all = TRUE)

# Frequency histogram
dt[!is.na(cel_freq)][, .(gene, crude_freq, cel_freq)] %>% 
  melt.data.table(id.vars = "gene") %>% 
  ggplot2::ggplot(aes(x = value, fill = variable)) + 
  ggplot2::geom_histogram(position = position_dodge()) + 
  ggplot2::geom_vline(xintercept = 0.95, lty = 2) +
  ggplot2::scale_fill_discrete(name = "Estimate type", labels = c("crude", "celebrimbor")) +
  ggplot2::labs(x = "Gene frequency (%)", y = "Number of genes") +
  cowplot::theme_minimal_grid() +
  ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1), labels = seq(0, 100, 10))
