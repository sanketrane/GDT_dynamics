functions{
  real theta_spline(real t,
    real nu,               // nu is the rate of decline of source influx assuming that its proportional to the counts of source compartment
    real theta0) {     // theta0 gives the initial counts of the source compartment

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

  real[] khm(real t, real[] y, real[] parms, real[] rdata, int[] idata) {
     real dydt[4];             // system of ODEs

     real psi      = parms[1];                          // daily input from source as a frcation of total pool.
     real f_fast   = parms[2];                          // fraction of fast cells within the source influx.
     real lambda_f = parms[3];                          // net rate of loss of fast subset.
     real lambda_s = parms[4];                          // net rate of loss of slow subset.

     real theta0 = rdata[1];
     real nu     = rdata[2];
     real chiEst = rdata[3];
     real qEst   = rdata[4];

     dydt[1] = (psi * theta_spline(t, nu, theta0) * Chi_spline(t, chiEst, qEst) * f_fast) - lambda_f * y[1];
     dydt[2] = (psi * theta_spline(t, nu, theta0) * Chi_spline(t, chiEst, qEst) * (1-f_fast)) - lambda_s * y[2];
     dydt[3] = (psi * theta_spline(t, nu, theta0) * (1 - Chi_spline(t, chiEst, qEst)) * f_fast) -  lambda_f  * y[3];
     dydt[4] = (psi * theta_spline(t, nu, theta0) * (1 - Chi_spline(t, chiEst, qEst)) * (1-f_fast)) - lambda_s * y[4];
     return dydt;
 }

 // functions for transformation of fractions in (0,a), where a >=1
 real[] asinsqrt_array(real[] x){
   real a = 1.2;

   real ans = Arcsin(a^(-0.5) * x^(0.5));
   return ans;
 }

 real asinsqrt_real(real x){
   real a = 1.2;

   real ans = asin(a^(-0.5) * x^(0.5));
   return ans;
 }

 real asinsqrt_inv(real x){
   real a = 1.2;

   real ans = a * (sin(x))^2;
   return ans;
 }
}

data{
  int<lower  = 1> numObs;                                // length of the vector of all the timepoints in a data set.
  int<lower  = 1> num_index;                             // length of the vector of the unique timepoints.
  real<lower = 0> solve_time[num_index];                 // unique time points to solve ODEs.
  int<lower  = 0> time_index[numObs];                    // indext to track timepoints between ODE solution and the observations.
  int dpBMT[numObs];                                     // days post BMT to feed in Chi_spline function for fd normalisation
  real Nd_0;                                             // Nd_0 is the count of host population at BMT i.e time zero => Nd_0 = 0.
  real<lower = 0> counts[numObs];
  real<lower = 0> Nfd[numObs];
  int<lower  = 1> numPred;                               // length of the vector used for predction timepoints
  real ts_pred[numPred];                                 // timpoints for predction
  real theta0;
  real nu;
  real chiEst;
  real qEst;
  }

transformed data{
  real y1[numObs];                  // data set 1
  real y2[numObs];                  // data set 2
  real rdata[4];                    // model input (fixed parameters) type real.
  int idata[0];                     // model input (fixed parameters) type integers.

  rdata[1] = theta0;                // theta0 is the count of source population at BMT i.e time zero.
  rdata[2] = nu;                    // nu is the rate of change source population with time.
  rdata[3] = chiEst;
  rdata[4] = qEst;

  y1 = log(counts);                 // transforming cell counts of donor compartments to feed in to ODEs
  y2 = asinsqrt_array(Nfd);                       // untransfored cell counts of donor fractions normalised to source chimerism to feed in to ODEs
}

parameters{
  real<lower = 0, upper=1> psi;
  real<lower = 0, upper=1> f_fast;
  real<lower = 0> lambda_f;
  real<lower = 0> lambda_s;

  real y0_Log;                           // log transformed y0
  real<lower = 0, upper=1> alpha;        // alpha is fraction of fast cells within the Target population.

  real<lower = 0> sigma1;                  // variance for data set 1
  real<lower = 0> sigma2;                  // variance for data set 2
  }

transformed parameters{
  real y_hat[num_index, 4];          // array assigned to ODE solutions (dims in the bracket).
  real y1_mean[numObs];
  real y2_mean[numObs];
  real parms[4];
  real init_cond[4];
  real y0;

  y0 = exp(y0_Log);

  init_cond[1] = Nd_0;               // Nd_0 is the count of host population at BMT i.e time zero => Nd_0 = 0.
  init_cond[2] = Nd_0 ;
  init_cond[3] = y0 * alpha;
  init_cond[4] = y0 * (1 - alpha);   // y0 is the count of host population at BMT i.e time zero.
  parms[1] = psi;
  parms[2] = f_fast;
  parms[3] = lambda_f;
  parms[4] = lambda_s;

  // Ode Solver
  y_hat[1, ] = init_cond;
  y_hat[2:num_index, ] = integrate_ode_rk45(khm, init_cond, solve_time[1], solve_time[2:num_index], parms, rdata, idata);

  for (i in 1:numObs){
    // total counts
    y1_mean[i] = y_hat[time_index[i], 1] + y_hat[time_index[i], 2] + y_hat[time_index[i], 3] + y_hat[time_index[i], 4];

    // donor fractions normalised with chimerism in the source
    y2_mean[i] = (y_hat[time_index[i], 1] + y_hat[time_index[i], 2])/(y1_mean[i] * Chi_spline(dpBMT[i], chiEst, qEst));
  }
}

model{
  psi ~ normal(0.5, 0.25);
  f_fast ~ normal(0.5, 0.25);
  lambda_f ~ normal(0.01, 1);
  lambda_s ~ normal(0.01, 1);

  y0_Log ~ normal(11, 2);
  alpha ~ normal(0.5, 0.25);

  sigma1 ~ cauchy(0.1, 2);
  sigma2 ~ cauchy(0.1, 2);

  y1 ~ normal(log(y1_mean), sigma1);
  y2 ~ normal(asinsqrt_array(y2_mean), sigma2);
}

generated quantities{
  real y_hat_pred[numPred, 4];
  real y1_mean_pred[numPred];
  real y2_mean_pred[numPred];
  real countspred[numPred];
  real fdpred[numPred];
  vector[numObs] log_lik;
  vector[numObs] log_lik1;
  vector[numObs] log_lik2;

  y_hat_pred[1, ] = init_cond;
  y_hat_pred[2:numPred, ] = integrate_ode_rk45(khm, init_cond, ts_pred[1], ts_pred[2:numPred], parms, rdata, idata);

  // Initial conditions for total counts, fd and fractions of ki67host and ki67 donor
  y1_mean_pred[1] = y0;
  countspred[1] = exp(normal_rng(log(y1_mean_pred[1]), sigma1));

  y2_mean_pred[1] = Nd_0;
  fdpred[1] = Nd_0;

  for (i in 2:numPred){
    y1_mean_pred[i] = y_hat_pred[i, 1] + y_hat_pred[i, 2] + y_hat_pred[i, 3] + y_hat_pred[i, 4];
    countspred[i] = exp(normal_rng(log(y1_mean_pred[i]), sigma1));

    y2_mean_pred[i] = (y_hat_pred[i, 1] + y_hat_pred[i, 2]) / (y1_mean_pred[i] * Chi_spline(ts_pred[i], chiEst, qEst));
    fdpred[i] = asinsqrt_inv(normal_rng(asinsqrt_real(y2_mean_pred[i]), sigma2));
  }

  // calculating the log predictive accuracy for each point
  for (n in 1:numObs) {
    log_lik1[n] = normal_lpdf(y1[n] | log(y1_mean[n]), sigma1);
    log_lik2[n] = normal_lpdf(y2[n] | asinsqrt_real(y2_mean[n]), sigma2);
    log_lik[n] = log_lik1[n] + log_lik2[n];
  }
}
