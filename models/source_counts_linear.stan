functions{
  real theta_spline(real t,
    real mu,               // nu is the rate of decline of source influx assuming that its proportional to the counts of source compartment
    real theta0) {     // theta0 gives the initial counts of the source compartment

      real theta;

      theta = theta0 * exp(-mu * t);
      return theta;
   }
 }


data{
  int<lower = 1> numObs;
  int<lower = 0> Time[numObs];
  int<lower = 1> numPred;
  int<lower = 0>  ts_pred[numPred];
  real<lower = 0> counts[numObs];
}

transformed data{
  real y[numObs];

  y = log(counts);
}

parameters{
  real y0Log;
  real mu;
  real<lower = 0> sigma;
}

transformed parameters{
  real ymean[numObs];
  real y0 = exp(y0Log);

  for(i in 1:numObs){
    ymean[i] = theta_spline(Time[i], mu, y0);
  }
}

model{
  mu ~ normal(0, 1);
  y0Log ~ normal(11, 2);
  sigma ~ cauchy(0.5, 2);

  y ~ normal(log(ymean), sigma);
}

generated quantities{
  real ymean_pred[numPred];
  real countspred[numPred];

  for(i in 1:numPred){
    ymean_pred[i] = (ts_pred[i] * mu) + y0;
    countspred[i] = exp(normal_rng(log(ymean_pred[i]), sigma));
  }
}
