data {
  int<lower=0> n;                  // length of each subcomponent -> total dim = 2n
  real<lower=0> s_hi;              // higher sd
  real<lower=0,upper=s_hi> s_lo;   // lower sd
}

parameters {
  vector[n] x_hi;                  // high sd component
  vector[n] x_lo;                  // low sd component
}

model {
  x_hi ~ normal(0, s_hi); 
  x_lo ~ normal(0, s_lo);
}

