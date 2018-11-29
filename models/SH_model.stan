functions{
  real theta_spline(real t,
    real nu,               // nu is the rate of decline of source influx assuming that its proportional to the counts of source compartment
    real theta0) {         // theta0 gives the initial counts of the source compartment

      real theta;

      theta = theta0 * exp(-nu * t);
      return theta;
   }

   real Chi_spline( real t,
     real chiEst,         // chiEst is the level if stabilised chimerism in the source compartment
     real qEst) {         // qEst is the rate with which cimerism chnages in the source compartment

       real chi;

       if (t < 0){
         chi = 0;                       // conditioning the function to adapt to the timepoints before BMT
       } else {
         chi = chiEst * (1 - exp(-qEst * t));
       }
       return chi;
    }

  real[] shm(real t, real[] y, real[] parms, real[] rdata, int[] idata) {
     real dydt[2];             // system of ODEs

     real psi    = parms[1];
     real lambda = parms[2];

     real theta0 = rdata[1];
     real nu     = rdata[2];
     real chiEst = rdata[3];
     real qEst   = rdata[4];

     dydt[1] = (psi * theta_spline(t, nu, theta0) * Chi_spline(t, chiEst, qEst)) - lambda * y[1];
     dydt[2] = (psi * theta_spline(t, nu, theta0) * (1 - Chi_spline(t, chiEst, qEst))) - lambda * y[2];
     return dydt;
 }

 // functions for transformation of fractions in (0,a), where a >=1
 real[] asinsqrt_array(real[] x){
   int ndims = size(x);
   real answer[ndims];
   real a = 1.2;

   for (i in 1: ndims){
     answer[i] = asin(sqrt(x[i])/sqrt(a));
   }
   return answer;
 }

 real asinsqrt_real(real x){
   real a = 1.2;

   real answer = asin(sqrt(x)/sqrt(a));
   return answer;
 }

 real asinsqrt_inv(real x){
   real a = 1.2;

   real answer = a * (sin(x))^2;
   return answer;
 }
}

data{
  int<lower  = 1> numObs;
  int<lower  = 1> num_index;
  real<lower = 0> solve_time[num_index];
  int<lower  = 0> time_index[numObs];
  int dpBMT[numObs];
  real Nd_0;
  real<lower = 0> counts[numObs];
  real<lower = 0> Nfd[numObs];
  int<lower  = 1> numPred;
  real ts_pred[numPred];
  real theta0;
  real nu;
  real chiEst;
  real qEst;
  }

transformed data{
  real y1[numObs];
  real y2[numObs];
  real rdata[4];
  int idata[0];

  rdata[1] = theta0;
  rdata[2] = nu;
  rdata[3] = chiEst;
  rdata[4] = qEst;

  y1 = log(counts);                 // transforming cell counts of donor compartments to feed in to ODEs
  y2 = (Nfd);                       // untransfored cell counts of donor fractions normalised to source chimerism to feed in to ODEs
}

parameters{
  real y0_Log;
  real<lower = 0, upper=1> psi;
  real<lower = 0> lambda;

  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  }

transformed parameters{
  real y_hat[num_index, 2];
  real y1_mean[numObs];
  real y2_mean[numObs];
  real parms[2];
  real init_cond[2];
  real y0;

  y0 = exp(y0_Log);

  init_cond[1] = Nd_0;
  init_cond[2] = y0 ;
  parms[1] = psi;
  parms[2] = lambda;

  y_hat[1, ] = init_cond;
  y_hat[2:num_index, ] = integrate_ode_rk45(shm, init_cond, solve_time[1], solve_time[2:num_index], parms, rdata, idata);

  for (i in 1:numObs){
    // total counts
    y1_mean[i] = y_hat[time_index[i], 1] + y_hat[time_index[i], 2];

    // donor fractions normalised with chimerism in the source
    y2_mean[i] = (y_hat[time_index[i], 1])/(y1_mean[i] * Chi_spline(dpBMT[i], chiEst, qEst));
  }
}

model{
  psi ~ normal(0.5, 0.25);
  y0_Log ~ normal(11, 2);
  lambda ~ normal(0.01, 1);

  sigma1 ~ cauchy(0.1, 1);
  sigma2 ~ cauchy(0.1, 1);

  y1 ~ normal(log(y1_mean), sigma1);
  y2 ~ normal((y2_mean), sigma2);
}

generated quantities{
  real y_hat_pred[numPred, 2];
  real y1_mean_pred[numPred];
  real y2_mean_pred[numPred];
  real countspred[numPred];
  real fdpred[numPred];
  vector[numObs] log_lik;
  vector[numObs] log_lik1;
  vector[numObs] log_lik2;

  y_hat_pred[1, ] = init_cond;
  y_hat_pred[2:numPred, ] = integrate_ode_rk45(shm, init_cond, ts_pred[1], ts_pred[2:numPred], parms, rdata, idata);

  // Initial conditions for total counts, fd and fractions of ki67host and ki67 donor
  y1_mean_pred[1] = y0;
  countspred[1] = exp(normal_rng(log(y1_mean_pred[1]), sigma1));

  y2_mean_pred[1] = Nd_0;
  fdpred[1] = Nd_0;

  for (i in 2:numPred){
    y1_mean_pred[i] = y_hat_pred[i, 1] + y_hat_pred[i, 2];
    countspred[i] = exp(normal_rng(log(y1_mean_pred[i]), sigma1));

    y2_mean_pred[i] = y_hat_pred[i, 1] / (y1_mean_pred[i] * Chi_spline(ts_pred[i], chiEst, qEst));
    fdpred[i] = (normal_rng((y2_mean_pred[i]), sigma2));
  }

  // calculating the log predictive accuracy for each point
  for (n in 1:numObs) {
    log_lik1[n] = normal_lpdf(y1[n] | log(y1_mean[n]), sigma1);
    log_lik2[n] = normal_lpdf(y2[n] | (y2_mean[n]), sigma2);
    log_lik[n] = log_lik1[n] + log_lik2[n];
  }
}
