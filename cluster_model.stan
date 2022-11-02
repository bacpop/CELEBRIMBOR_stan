
data {
  int <lower = 0> N;
  int <lower = 0> K;
  int <lower = 0> C;
  array[N, K] int dat;
  vector[N] comp;
  array[N] int clust;
}

parameters {
  matrix <lower = 0, upper = 1> [C, K] true_freq;
}

model {
  for(i in 1:N){
    dat[i, ] ~ bernoulli(true_freq[clust[i], ] * comp[i]);
  }
}
