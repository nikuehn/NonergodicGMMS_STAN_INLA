# ------------------------------------------------------------
# script to run stan models from R using cmdstanR
#
# see https://mc-stan.org/cmdstanr/ for more information
#
# note that for the paper, the models were run directly
# from the command line using cmdstan
# -----------------------------------------------------------
#clear variables
rm(list = ls())
# load packages

# function definitions
local.plot.field = function(field, mesh, xlim=c(0,10), ylim=c(0,10), ...){
  stopifnot(length(field) == mesh$n)
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  field.proj = inla.mesh.project(proj, field)
  n.col = 20
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, col = plasma(n.col), nlevel=n.col+1)
}

library(cmdstanr)
library(posterior)
library(bayesplot)

set_cmdstan_path('/Users/nico/GROUNDMOTION/SOFTWARE/cmdstan-2.26.1')
cmdstan_path()
cmdstan_version()

setwd('/Users/nico/GROUNDMOTION/PROJECTS/NONERGODIC/COMPARISON_REGRESSION/Git/NonergodicGMMS_STAN_INLA/')

# read in Stan data file, which contains all relevant data
period = 0.02
dat_stan <- rstan::read_rdump(sprintf('DATA/data_vcm_cellattn_nngp_ASK_%s.Rdata',period))

# define Stan model for NNGP
file <- file.path('STAN_CODE', 'gmm_NNGPeq_NNGPstat_cellattn.stan')
mod <- cmdstan_model(file)

# specify output directory to save csv files
# run only 50/50  samples and 2 chains here for speed
fit <- mod$sample(
  data = 'DATA/data_vcm_cellattn_nngp_ASK_0.02.Rdata',
  seed = 5618,
  chains = 2,
  iter_sampling = 50,
  iter_warmup = 50,
  refresh = 10,
  max_treedepth = 15,
  adapt_delta = 0.9,
  parallel_chains = 2,
  output_dir = 'SAMPLES_STAN/'
)
# check that fit is ok (will not be for 50 samples)
fit$cmdstan_diagnose()
fit$save_object(file = file.path('SAMPLES_STAN', 'fit_nngp_ask.RDS'))
fit <- readRDS(file.path('SAMPLES_STAN', 'fit_nngp_ask.RDS'))

# access samples
draws <- fit$draws()

summarise_draws(subset(draws, variable=c( "phi_0", "tau_0", "phi_S2S")))
summarise_draws(subset(draws, variable=c( "rho_", "theta_"), regex = TRUE))

# plot posterior distributions of aleatory parameters
mcmc_areas(
  draws,
  pars = c("phi_0", "tau_0", "phi_S2S"),
  prob = 0.68, # 80% intervals
  prob_outer = 0.9, # 99%
  point_est = "mean"
)


# run Stan model to predict even/station terms
load('RESULTS_INLA/mesh.Rdata')

# the following model reads in the posterior samples of
# spatially varying event terms, length scale, and standard deviation
# calculates predictive mean/standard deviation at new locations
# samples from this distribution, and calculates mean/sd at each location
file <- file.path('STAN_CODE', 'prediction_gp_posterior_rng.stan')
mod_pred <- cmdstan_model(file)

params <- as_draws_df(subset(draws, variable = c('rho_eq','theta_eq')))
# spatial varying event terms are in order of nearest neighbours
feq <- as_draws_df(subset(draws, variable = c('f_eq')))[,dat_stan$eq_nn2]

mesh_nodes <- mesh1$loc

data_list <- list(N = dat_stan$NEQ,
                  N_samples = length(params$rho_eq),
                  N_star = length(mesh_nodes[,1]),
                  rho = params$rho_eq,
                  theta = params$theta_eq,
                  Y = feq,
                  X = dat_stan$X_e,
                  X_star = mesh_nodes[,c(1,2)]
)

fit_pred <- mod_pred$sample(
  data = data_list,
  fixed_param = TRUE,
  iter_sampling = 1,
  iter_warmup = 0,
)
draws_pred <- fit_pred$draws()

mu_star <- as_draws_df(subset(draws_pred, variable = "mu_star"))
sd_star <- as_draws_df(subset(draws_pred, variable = "sd_star"))


### calculate predictions with point estimates
file <- file.path('STAN_CODE', 'prediction_gp_psi.stan')
mod_pred2 <- cmdstan_model(file)

cov_feq <- cov(feq)

data_list <- list(N = dat_stan$NEQ,
                  N_star = length(mesh_nodes[,1]),
                  rho = median(params$rho_eq),
                  theta = median(params$theta_eq),
                  Y = colMeans(feq),
                  psi = cov_feq,
                  X = dat_stan$X_e,
                  X_star = mesh_nodes[,c(1,2)]
)

fit_pred2 <- mod_pred2$sample(
  data = data_list,
  fixed_param = TRUE,
  iter_sampling = 1,
  iter_warmup = 0,
)
draws_pred2 <- fit_pred2$draws()

tmp <- subset_draws(draws_pred2, "mu_star")[1:mesh1$n]
local.plot.field(tmp[1:mesh1$n], mesh1, xlim = c(-126,-114), ylim = c(30,42))

tmp <- subset_draws(draws_pred2, "sd_star")[1:mesh1$n]
local.plot.field(tmp[1:mesh1$n], mesh1, xlim = c(-126,-114), ylim = c(30,42))
