functions {
  matrix updatePrior(int[] pos, real a0, real m) {
    
    matrix[num_elements(pos), 10] mu = rep_matrix(0, num_elements(pos), 10);
    row_vector[10] a = rep_row_vector(a0, 10); 
    
    mu[1,1:10] = a/sum(a);
    for(i in 2:num_elements(pos)) {
      if(pos[i-1] > 0) {
        a[pos[i-1]] = a[pos[i-1]] + 1;
      }
      a = a*m + (1-m)*rep_row_vector(a0, 10); 
      mu[i] = a/sum(a);
    }
    return mu;
  } 
}

data {
  int N; //n number of trials
  real rt[N]; // response times
  int dpos[N]; // distractor position
  int tpos[N]; // target position
}

parameters {
  real<lower=0> a0;
  real<lower=0,upper=1> mem;
  real<lower=0> k;
  real t_a;
  real k_d;
  real d_a;
  real sigma;
}

transformed parameters {
  
  // Response time predictions
  matrix[N,10] mu = updatePrior(dpos, a0, mem);
  vector[N] t_prob;
  vector[N] d_prob  = rep_vector(0,N);
  vector[N] dist_present = rep_vector(0,N);
  vector[N] pred_rt;
  
  for(i in 1:N) {
    t_prob[i] = mu[i, tpos[i]];
    if (dpos[i] > 0) {
      d_prob[i] = mu[i, dpos[i]];
      dist_present[i] = 1;
    }
  }
  
  pred_rt = k + t_a*t_prob + k_d*dist_present - d_a*d_prob; 
}

model {
  // Priors
  mem ~ uniform(0,1); 
  a0 ~ normal(1, 30);
  k ~ normal(600, 300);
  t_a ~ normal(400, 200);
  k_d ~ normal(200, 100);
  d_a ~ normal(700, 300);

  // Likelihood
  for(i in 1:N) {
    rt[i] ~ normal(pred_rt[i], sigma);
  }
}
