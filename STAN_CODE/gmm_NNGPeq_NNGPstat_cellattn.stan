/*********************************************
Stan program to obtain VCM parameters
lower dimensions is used (event terms/station terms)

This model explicitly estimates the latent (uncorrelated) event terms and station terms

event constant is correlated with Matern 1/2
station constant is correlated with Matern 1/2

spatially correlated event and station terms are modeled as NNGPs

likelihood of NGP is taken from https://mc-stan.org/users/documentation/case-studies/nngp.html

anelastic attenuaion is modeled with a cell-specific attenuation
 ********************************************/

functions{
  real nngp_w_lpdf(vector w, real sigmasq, real phi, matrix NN_dist,
                       matrix NN_distM, int[,] NN_ind, int N, int M){

      vector[N] V;
      vector[N] I_Aw = w;
      int dim;
      int h;

      for (i in 2:N) {

          matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
          iNNdistM;
          matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
          iNNCholL;
          vector[ i < (M + 1)? (i - 1) : M] iNNcorr;
          vector[ i < (M + 1)? (i - 1) : M] v;
          row_vector[i < (M + 1)? (i - 1) : M] v2;

          dim = (i < (M + 1))? (i - 1) : M;

          if(dim == 1){iNNdistM[1, 1] = 1;}
          else{
              h = 0;
              for (j in 1:(dim - 1)){
                  for (k in (j + 1):dim){
                      h = h + 1;
                      iNNdistM[j, k] = exp(- phi * NN_distM[(i - 1), h]);
                      iNNdistM[k, j] = iNNdistM[j, k];
                  }
              }
              for(j in 1:dim){
                  iNNdistM[j, j] = 1 + 1e-9;
              }
          }

          iNNCholL = cholesky_decompose(iNNdistM);
          iNNcorr = to_vector(exp(- phi * NN_dist[(i - 1), 1:dim]));

          v = mdivide_left_tri_low(iNNCholL, iNNcorr);

          V[i] = 1 - dot_self(v);

          v2 = mdivide_right_tri_low(v', iNNCholL);

          I_Aw[i] = I_Aw[i] - v2 * w[NN_ind[(i - 1), 1:dim]];

      }
      V[1] = 1;
      return - 0.5 * ( 1. / sigmasq * dot_product(I_Aw, (I_Aw ./ V)) +
                      sum(log(V)) + N * log(sigmasq));
  }
}

data {
  int N;  // number of records
  int NEQ;  // number of earthquakes
  int NSTAT;  // number of stations
  int NCELL;  // number of cells
  int<lower=1> NN; // number of nearest neghbours

  matrix[N,NCELL] RC; // for each record, distance within each cell

  vector[N] mu_rec; // median predictions for each record

  vector[N] Y; // ln ground-motion value

  // parameters for NGPs
  int NN_ind_eq[NEQ - 1, NN];
  matrix[NEQ - 1, NN] NN_dist_eq;
  matrix[NEQ - 1, (NN * (NN - 1) / 2)] NN_distM_eq;

  int NN_ind_stat[NSTAT - 1, NN];
  matrix[NSTAT - 1, NN] NN_dist_stat;
  matrix[NSTAT - 1, (NN * (NN - 1) / 2)] NN_distM_stat;

  int<lower=1,upper=NEQ> eq[N]; // event id (in numerical order from 1 to last)
  int<lower=1,upper=NSTAT> stat[N]; // station id (in numerical order from 1 to last)

  int<lower=1,upper=NEQ> eq_nn2[NEQ];  // event id for nearest neighbours
  int<lower=1,upper=NSTAT> stat_nn2[NSTAT]; // saion id for nearest neghbour
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

  vector[NEQ] eqterm;
  vector[NSTAT] statterm;

  vector[NEQ] f_eq;
  vector[NSTAT] f_stat;

  vector<lower=0>[NCELL] c_ca; // cell-specific attenuation
}

transformed parameters {
  real lambda_stat = inv(rho_stat);
  real lambda_eq = inv(rho_eq);
}

model {
  vector[N] mu_rec2;

  intercept ~ normal(0,0.1);

  phi_0 ~ lognormal(ln05, sig);
  tau_0 ~ lognormal(ln05, sig);
  phi_S2S ~ lognormal(ln05, sig);

  rho_eq ~ inv_gamma(3,0.5);
  rho_stat ~ inv_gamma(3,0.5);
  theta_eq ~ exponential(20);
  theta_stat ~ exponential(20);

  eqterm ~ normal(0,tau_0);
  statterm ~ normal(0,phi_S2S);

  // priors anelastic attenuation
  sigma_ca ~ exponential(100);
  mu_ca ~ lognormal(-4.6,0.5);  // prior for global attenuation parameter
  c_ca ~ normal(mu_ca,sigma_ca);  // prior for regional attenuatio

  // latent variable NNGP terms
  f_eq ~ nngp_w(square(theta_eq), lambda_eq, NN_dist_eq, NN_distM_eq, NN_ind_eq, NEQ, NN);
  f_stat ~ nngp_w(square(theta_stat), lambda_stat, NN_dist_stat, NN_distM_stat, NN_ind_stat, NSTAT, NN);

  mu_rec2 = mu_rec + intercept - RC * c_ca + f_eq[eq_nn2[eq]] + f_stat[stat_nn2[stat]] + eqterm[eq] + statterm[stat];

  Y ~ normal(mu_rec2,phi_0);

}
