/* *****************************************************************
Program to calculate mean and sandard deviation of predictions at new locatios

prediction is based on point estimates
takes uncerainty of latent function values into account via Psi matrix
* *****************************************************************/

data {
  int N;
  int N_star;

  real<lower=0> rho;
  real<lower=0> theta;
  vector[N] Y;
  matrix[N, N] psi;

  vector[2] X[N];
  vector[2] X_star[N_star];
}

transformed data {
  real delta = 1e-9;
}

parameters {

}

model {

}

generated quantities {
  vector[N_star] mu_star;
  vector[N_star] sd_star;

  {
    matrix[N, N] cov;
    // calculate covariance matrix
    for(i in 1:N) {
      for(j in (i + 1):N) {
        cov[i,j] = square(theta) * exp( - distance(X[i], X[j]) /rho);
        cov[j,i] = cov[i,j];
      }
      cov[i,i] = square(theta) + delta;
    }

    for(i in 1:N_star) {
      row_vector[N] kvec;
      for(j in 1:N) {
        kvec[j] = square(theta) * exp( - distance(X_star[i], X[j]) /rho);
      }
      row_vector[N] kvec_inv_cov = kvec / cov;
      mu_star[i] = (kvec / cov) * Y;
      sd_star[i] = sqrt(square(theta) - kvec_inv_cov * to_vector(kvec) + quad_form(psi, to_vector(kvec_inv_cov)));
    }
  }
}
