functions {

  /**
    * Return the Poisson-binomial log probability mass for the specified
  * count y and vector of probabilities theta.  The result is the log
  * probability of y successes in N = size(theta) independent
  * trials with probabilities of success theta[1], ..., theta[N].
  *
    * See:  https://en.wikipedia.org/wiki/Poisson_binomial_distribution
  *
    * @param y number of successes
  * @param theta vector of probabilities
  * @return Poisson-binomial log probability mass
  */
    real poisson_binomial_lpmf(int y, vector theta) {
      int N = rows(theta);
      matrix[N + 1, N + 1] alpha;
      vector[N] log_theta = log(theta);
      vector[N] log1m_theta = log1m(theta);
      
      if (y < 0 || y > N)
        reject("poisson binomial variate y out of range, y = ", y,
               " N = ", N);
      for (n in 1:N)
        if (theta[n] < 0 || theta[n] > 1)
          reject("poisson binomial parameter out of range,",
                 " theta[", n, "] =", theta[n]);
      
      if (N == 0)
        return y == 0 ? 0 : negative_infinity();
      
      // dynamic programming with cells
      // alpha[n + 1, tot + 1] = log prob of tot successes in first n trials
      alpha[1, 1] = 0;
      for (n in 1:N) {
        // tot = 0
        alpha[n + 1, 1] = alpha[n, 1] + log1m_theta[n];
        
        // 0 < tot < n
        for (tot in 1:min(y, n - 1))
          alpha[n + 1, tot + 1]
        = log_sum_exp(alpha[n, tot] + log_theta[n],
                      alpha[n, tot  + 1] + log1m_theta[n]);
        
        // tot = n
        if (n > y) continue;
        alpha[n + 1, n + 1] = alpha[n, n] + log_theta[n];
      }
      return alpha[N + 1, y + 1];
    }
  
  //   real partial_sum(array[] int y_slice,
  //                  int start, int end,
  //                  vector pars,
  //                  vector completeness) {
  //   real ll = 0;
  //   int n = num_elements(y_slice);
  //   
  //   for(i in 1:n){
  //     ll += poisson_binomial_lpmf(y_slice[i] | pars[start - 1 + i] * completeness);
  //   }
  //   
  //   return ll;
  // }
}

data {
  // Number of genome samples
  int <lower = 0> N;
  // Array of # of gene occurences
  int y;
  // Vector of genome sample completeness
  vector[N] completeness;
  // Beta prior parameters
  real a;
  real b;
}

parameters {
  real <lower = 0, upper = 1> true_freq;
}

model {
  vector[N] p = true_freq * completeness;
  true_freq ~ beta(a, b);
  target += poisson_binomial_lpmf(y | p);
}
