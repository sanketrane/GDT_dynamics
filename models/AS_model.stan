functions{
  real lambda(real a, real t, real lambda_0, real r_l, real init) {
    real ans;

    if (init >= 1){
      ans = (lambda_0 /r_l) * (exp(-r_l * (a - t)) - exp(-r_l * a));
    } else {
      ans = (lambda_0 /r_l) * (1 - exp(-r_l * a));
    }
    return ans;
  }

  real theta(real t, real N0, real p) {
    real ans;
    int ta = 40;
    real nu = 0.0014307;
    real theta0 = N0 * p/(exp(p * ta) -1);

    if (t >= 0){
      ans = theta0 * exp(-nu * t);
    } else {
      ans = 0;
    }
    return ans;
  }

  real Chi_spline( real t, real chiEst, real qEst) {
    // chiEst is the level if stabilised chimerism in the source compartment
    // qEst is the rate with which cimerism chnages in the source compartment
    real chi;

    if (t < 0){
      chi = 0;                       // conditioning the function to adapt to the timepoints before BMT
    } else {
      chi = chiEst * (1 - exp(-qEst * t));
    }
    return chi;
  }

  real g_age(real a, real N0, real p) {
    real ans;
    int ta = 40;
    real theta0 = N0 * p/(exp(p * ta) -1);

    if (a >= 0){
      ans = theta0 * exp(p * a);
    } else {
      ans = 0;
    }
    return ans;
  }


  real[] ode_func(real a,  real[] y, real[] parms, real[] rdata, int[] idata) {
   real lambda_0 = parms[1];
   real r_l = parms[2];
   real N0 = parms[3];
   real p_age = parms[4];

   real theta0 = rdata[1];           // initial value of theta function
   real nu = rdata[2];               // rate of change of theta with t
   real chiEst = rdata[3];
   real qEst = rdata[4];
   real dydt[3];

   int tau = idata[1];

   dydt[1] = theta(tau - a, N0, p_age) * Chi_spline(tau, chiEst,  qEst) * exp(- lambda(a, tau, lambda_0, r_l, 0));
   dydt[2] = theta(tau - a, N0, p_age) * (1 - Chi_spline(tau, chiEst,  qEst)) * exp(- lambda(a, tau, lambda_0, r_l, 0));
   dydt[3] = g_age(a - tau, N0, p_age) * exp(- lambda(a, tau, lambda_0, r_l, 1));

   return dydt;
  }

 real[] foreach_ode(real ts, real t0, real[] init_cond, real[] parms, real[] rdata, int x_i) {
   return to_array_1d(integrate_ode_rk45(ode_func, init_cond, t0, rep_array(ts, 1), parms, rdata, rep_array(x_i, 1)));
 }

 real[,] solve_ode(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] d_var, int num_index) {

   real y_hat[num_index, 3];

   y_hat[1] = init_cond;
   for (i in 2:num_index){
     y_hat[i] = foreach_ode(solve_time[i], solve_time[1], init_cond, parms, rdata, d_var[i]);
   }

  return y_hat;
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
  int tau_time[num_index];
  int dpBMT[numObs];
  real Nd_0;
  real<lower = 0> counts[numObs];
  real<lower = 0> Nfd[numObs];
  int<lower  = 1> numPred;
  int tau_pred[numPred];
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
  int idata[num_index];

  for (i in 1:num_index){
    idata[i] = tau_time[i];
  }

  rdata[1] = theta0;
  rdata[2] = nu;
  rdata[3] = chiEst;
  rdata[4] = qEst;

  y1 = log(counts);                 // transforming cell counts of donor compartments to feed in to ODEs
  y2 = asinsqrt_array(Nfd);                       // untransfored cell counts of donor fractions normalised to source chimerism to feed in to ODEs
}

parameters{
  real<lower = 0> lambda0;
  real r_l;
  real<lower = 0> y0_Log;
  real p_age;

  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  }

transformed parameters{
  real y_hat[num_index, 3];
  real y1_mean[numObs];
  real y2_mean[numObs];
  real y3_mean[numObs];
  real parms[4];
  real init_cond[3];
  real y0;

  y0 = exp(y0_Log);

  init_cond[1] = Nd_0;
  init_cond[2] = Nd_0 ;
  init_cond[3] = y0 ;
  parms[1] = lambda0;
  parms[2] = r_l;
  parms[3] = y0;
  parms[4] = p_age;


  y_hat = solve_ode(solve_time, init_cond, parms, rdata, idata, num_index);

  for (i in 1:numObs){

    // total counts
    y1_mean[i] = y_hat[time_index[i], 1] + y_hat[time_index[i], 2] + y_hat[time_index[i], 3];

    // donor fractions normalised with chimerism in the source
    y2_mean[i] = (y_hat[time_index[i], 1])/(y1_mean[i] * Chi_spline(dpBMT[i], chiEst, qEst));
  }
}

model{
  lambda0 ~ normal(0.01, 1);
  r_l ~ normal(0.001, 1);
  y0_Log ~ normal(11, 2);
  p_age ~ normal(0.001, 1);

  sigma1 ~ cauchy(0.1, 1);
  sigma2 ~ cauchy(0.1, 1);

  y1 ~ normal(log(y1_mean), sigma1);
  y2 ~ normal(asinsqrt_array(y2_mean), sigma2);
}

generated quantities{
  real y_hat_pred[numPred, 3];
  real y1_mean_pred[numPred];
  real y2_mean_pred[numPred];
  real countspred[numPred];
  real fdpred[numPred];
  vector[numObs] log_lik;
  vector[numObs] log_lik1;
  vector[numObs] log_lik2;

  y_hat_pred = solve_ode(ts_pred, init_cond, parms, rdata, tau_pred, numPred);

  // Initial conditions for total counts, fd and fractions of ki67host and ki67 donor
  y1_mean_pred[1] = y0;
  countspred[1] = exp(normal_rng(log(y1_mean_pred[1]), sigma1));

  y2_mean_pred[1] = Nd_0;
  fdpred[1] = Nd_0;

  for (i in 2:numPred){
    y1_mean_pred[i] = y_hat_pred[i, 1] + y_hat_pred[i, 2] + y_hat_pred[i, 3];
    countspred[i] = exp(normal_rng(log(y1_mean_pred[i]), sigma1));

    y2_mean_pred[i] = y_hat_pred[i, 1] / (y1_mean_pred[i] * Chi_spline(tau_pred[i], chiEst, qEst));
    fdpred[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y2_mean_pred[i]), sigma2));
  }

  // calculating the log predictive accuracy for each point
  for (n in 1:numObs) {
    log_lik1[n] = normal_lpdf(y1[n] | log(y1_mean[n]), sigma1);
    log_lik2[n] = normal_lpdf(y2[n] | asinsqrt_real(y2_mean[n]), sigma2);
    log_lik[n] = log_lik1[n] + log_lik2[n];
  }
}
