functions {
    // compute the mean, which is the analytic solution of the ODE
    // use the identity
    //   (exp(-bt) - exp(-dt))/(d-b) = (exp(-mt) - exp(-Mt))/D // M=max{d,b}, m=min{d,b}, D=M-m
    //     = exp(-mt)(1 - exp(-Dt))/D
    //     = -exp(-mt)expm1(-Dt)/D                             // mt>0 and Dt > 0 => avoids exp(>>>0)
    // Also
    //   expm1(-Dt)/D ~ (-Dt + D^2t^2/2)/D
    //     =t(Dt/2-1)
    real get_mu(real tmt0, real km0, real beta, real delta){
        if (tmt0 <= 0.0){
            return 0.0;
        }
        real m,D;
        if (delta>beta){
            m = beta;
            D = delta-beta;
        } else {
            m = delta;
            D = beta-delta;
        }
        real C = D < machine_precision() ? tmt0*(D*tmt0/2-1) : -expm1(-D*tmt0)/D;
        return km0 * exp(-tmt0*m) * C;
    }
}
data {
    int <lower=0> N; // number of observations
    array[N] real<lower=0> ts; // time of the observation
    array[N] real ys; // observed value
}
parameters {
    real<lower=-2,upper=1> lt0;
    real<lower=-5,upper=5> lkm0;
    real<lower=-5,upper=5> lbeta;
    real<lower=-5,upper=5> ldelta;
    real<lower=-2,upper=2> lsigma;
}
transformed parameters{
    // real t0, km0, beta, delta, sigma;
    real t0    = pow(10, lt0);
    real km0   = pow(10, lkm0);
    real beta  = pow(10, lbeta);
    real delta = pow(10, ldelta);
    real sigma = pow(10, lsigma);
}
model {
    // Priors are all uniform so they are implicit
    // Likelihood:
    //     y_i|params ~indep N(mu_i, sigma)
    // with
    //     mu_i = get_mu(tmt0, km0, beta, delta, sigma)
    for (i in 1:N) {
        ys[i] ~ normal(get_mu(ts[i] - t0, km0, beta, delta), sigma);
    }
}