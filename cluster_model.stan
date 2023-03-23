functions {
  
  // real partial_sum(array[] int y_slice,
  //                  int start, int end,
  //                  vector completion,
  //                  vector freqs) {
  //   return bernoulli_logit_lpmf(y_slice | completion[start:end] .* freqs);
  // }
  
  real partial_sum(array[] int y_slice,
                   int start, int end,
                   vector pars,
                   vector completion) {
    real completion = 0.95;
    int NN = 1571;
    return binomial_lpmf(y_slice | NN, completion * freqs[start:end]);
  }
}

data {
  // Number of genomes
  int <lower = 0> N;
  // Number of genes
  int <lower = 0> K;
  // int <lower = 0> C;
  // array[N, K] int dat;
  array[K] int y;
  // array[N * K] int genome;
  // array[N * K] int gene;
  // vector[1] completion;
  real completion;
  vector[N] comp;
  // array[N] int clust;
}

parameters {
  // matrix <lower = 0, upper = 1> [C, K] true_freq;
  vector <lower = 0, upper = 1> [K]  true_freq;
}

// transformed parameters {
//   vector [N * K] beta  = true_freq;
// }

model {
  target += reduce_sum(partial_sum, y, 1, comp, true_freq);
  // for(i in 1:N){
    // dat[i, ] ~ bernoulli(true_freq[clust[i], ] * comp[i]);
    // dat[i, ] ~ bernoulli(true_freq * comp[i]);
  // }
  // true_freq ~ beta(0.5, 0.5);
}
