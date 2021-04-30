/* *****************************************************************
Program to calculate mean and sandard deviation of predictions at ew locatios
based o posterior samples

for each posterior sample, predicive mean/sd are calculated
from this distribution a random sample is drawn, and then mean/sd are calculated
* *****************************************************************/

data {
  int N_samples;
  int N;
  int N_star;

  vector[N_samples] rho;
  vector[N_samples] theta;
  vector[N] Y[N_samples];

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
    vector[N_star] mu_rng[N_samples];
    for(k in 1:N_samples) {
      matrix[N, N] cov;
      // calculate covariance matrix
      for(i in 1:N) {
        for(j in (i + 1):N) {
          cov[i,j] = square(theta[k]) * exp( - distance(X[i], X[j]) /rho[k]);
          cov[j,i] = cov[i,j];
        }
        cov[i,i] = square(theta[k]) + delta;
      }

      for(i in 1:N_star) {
        row_vector[N] kvec;
        real mu;
        real sigma;
        for(j in 1:N) {
          kvec[j] = square(theta[k]) * exp( - distance(X_star[i], X[j]) /rho[k]);
        }
        mu = (kvec / cov) * Y[k];
        sigma = sqrt(square(theta[k]) - (kvec / cov) * to_vector(kvec));
        mu_rng[k,i] = normal_rng(mu, sigma);
      }
    }

    for(i in 1:N_star) {
      mu_star[i] = mean(mu_rng[:,i]);
      sd_star[i] = sd(mu_rng[:,i]);
    }
  }

}
