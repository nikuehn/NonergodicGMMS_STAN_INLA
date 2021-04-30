/*********************************************
Stan program to obtain VCM parameters
lower dimensions is used (event terms/station terms)

This model explicitly estimates the latent (uncorrelated) event terms and station terms

event constant is correlated with Matern 1/2
station constant is correlated with Matern 1/2

anelastic attenuaion is modeled with a cell-specific attenuation
 ********************************************/

data {
  int N;  // number of records
  int NEQ;  // number of earthquakes
  int NSTAT;  // number of stations
  int NCELL;  // number of cells

  matrix[N,NCELL] RC; // for each record, distance within each cell

  vector[N] mu_rec;         // median predictions for each record

  vector[N] Y; // ln ground-motion value

  int<lower=1,upper=NEQ> eq[N]; // event id (in numerical order from 1 to last)
  int<lower=1,upper=NSTAT> stat[N]; // station id (in numerical order from 1 to last)

  vector[2] X_e[NEQ];  // event coordinate for each record
  vector[2] X_s[NSTAT];  // station coordinate for each record

}

transformed data {
  real delta = 1e-9;
  real ln05 = log(0.5);
  real sig = 0.5;
}

parameters {
  real intercept;

  real<lower=0> phi_0;  // phi_0 - remaining aleatory variability of within-event residuals
  real<lower=0> tau_0;  // tau_0 - remaining aleatory variability of between-event residuals
  real<lower=0> phi_S2S;  // phi_s2s - remaining aleatory variability of between-station residuals

  real<lower=0> rho_eq; // length scale for spatially correlated event terms
  real<lower=0> theta_eq; // std of spatially correlated event terms
  real<lower=0> rho_stat; // length scale for spatially correlated station terms
  real<lower=0> theta_stat; // std of spatially correlated station terms

  real<lower=0> mu_ca;       // global (mean) attenuation
  real<lower=0> sigma_ca;   // std of cell-specific attenuation

  vector[NEQ] z_eq;
  vector[NSTAT] z_stat;

  vector[NEQ] eqterm;
  vector[NSTAT] statterm;

  vector<lower=0>[NCELL] c_ca; // cell-specific attenuation
}

transformed parameters {
  vector[NEQ] f_eq;
  vector[NSTAT] f_stat;

  // latent variable event contributions to GP
  {
    matrix[NEQ,NEQ] cov_eq;
    matrix[NEQ,NEQ] L_eq;

    for(i in 1:NEQ) {
      for(j in i:NEQ) {
        real d_e;
        real c_eq;
  
        d_e = distance(X_e[i],X_e[j]);
  
        c_eq = (theta_eq^2 * exp(-d_e/rho_eq));
  
        cov_eq[i,j] = c_eq;
        cov_eq[j,i] = cov_eq[i,j];
      }
      cov_eq[i,i] = cov_eq[i,i] + delta;
    }

    L_eq = cholesky_decompose(cov_eq);
    f_eq = L_eq * z_eq;
  }

  // latent variable station contributions to GP
  {
    matrix[NSTAT,NSTAT] cov_stat;
    matrix[NSTAT,NSTAT] L_stat;

    for(i in 1:NSTAT) {
      for(j in i:NSTAT) {
        real d_s;
        real c_stat;
  
        d_s = distance(X_s[i],X_s[j]);
  
        c_stat = (theta_stat^2 * exp(-d_s/rho_stat));
  
        cov_stat[i,j] = c_stat;
        cov_stat[j,i] = cov_stat[i,j];
      }
      cov_stat[i,i] = cov_stat[i,i] + delta;
    }

    L_stat = cholesky_decompose(cov_stat);
    f_stat = L_stat * z_stat;
  }
}

model {
  vector[N] mu_rec2;
  vector[N] mu_rec3;

  intercept ~ normal(0,0.1);

  phi_0 ~ lognormal(ln05, sig); // based on 20% reduction compared to ergodic
  tau_0 ~ lognormal(ln05, sig); // based on 20% reduction compared to ergodic
  phi_S2S ~ lognormal(ln05, sig); // based on 20% reduction compared to ergodic

  rho_eq ~ inv_gamma(3,0.5);
  rho_stat ~ inv_gamma(3,0.5);
  theta_eq ~ exponential(20);
  theta_stat ~ exponential(20);

  z_eq ~ std_normal();
  z_stat ~ std_normal();

  eqterm ~ normal(0, tau_0);
  statterm ~ normal(0, phi_S2S);

  // priors anelastic attenuation
  sigma_ca ~ exponential(100);
  mu_ca ~ lognormal(-4.6,0.5);  // prior for global attenuation parameter
  c_ca ~ normal(mu_ca,sigma_ca);  // prior for regional attenuation

  mu_rec2 = mu_rec + intercept - RC * c_ca + f_eq[eq] + f_stat[stat] + eqterm[eq] + statterm[stat];
  Y ~ normal(mu_rec2,phi_0);

}
