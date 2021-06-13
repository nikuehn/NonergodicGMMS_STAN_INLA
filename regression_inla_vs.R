# ------------------------------------------------------------
# script to calculate a Bayesian VCM on NGA West 2 using INLA
#
# for information on R-INLA (installation and so on), https://www.r-inla.org/
#
# this notebook performs the regression, and shws how one can
# assess results (get a summary, get posterior distributions)
# -----------------------------------------------------------
#clear variables
rm(list = ls())
# load packages

library(INLA)
library(fields)
library(viridisLite)

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

slcpo <- function(m, na.rm = TRUE) {
  - sum(log(m$cpo$cpo), na.rm = na.rm)
}

setwd('/Users/nico/GROUNDMOTION/PROJECTS/NONERGODIC/COMPARISON_REGRESSION/Git/NonergodicGMMS_STAN_INLA/')

# -----------------------------------------------------------
# Data

period = '0.02'

################################################################
############### ASK ############################################
################################################################
# read in Stan data file, which contains all relevant data
dat_stan <- rstan::read_rdump(sprintf('DATA/data_vcm_cellattn_nngp_ASK_%s.Rdata',period))

coords_stat <- dat_stan$X_s
coords_eq <- dat_stan$X_e

n_rec <- dat_stan$N
n_eq <- dat_stan$NEQ
n_stat <- dat_stan$NSTAT
n_cell <- dat_stan$NCELL
cell_idx <- dat_stan$cell_idx

eqid <- dat_stan$eq
statid <- dat_stan$stat
mu_rec <- dat_stan$mu_rec # 
Y <- dat_stan$Y
resid <- Y - mu_rec
RC <- dat_stan$RC
R <- rowSums(RC)
lnVS <- log(dat_stan$VS/400)


# -------------------------------------------
# --------- linear model --------------------
# estimate ergodic linear model for comparison
# model contains a constant, and even/station terms
# -------------------------------------------

df.covar <- data.frame(intercept = 1,
                       R = R,
                       eq = eqid,
                       stat = statid,
                       resid = resid)
df.covar$idx_cell <- 1:nrow(df.covar)

# Prior on the fixed effects
prior.fixed <- list(mean.intercept = 0, prec.intercept = 100,
                    mean = (list(R=-0.01, default=0)), prec = (list(R=10000, default=1)))

# prior on tau and phi_S2S, phi_SS
prior_prec_tau <- list(prec = list(prior = "loggamma", param = c(2, 0.5)))
prior_prec_phiS2S <- list(prec = list(prior = "loggamma", param = c(2, 0.5)))
prior_prec <- list(prec = list(prior = "loggamma", param = c(2, 0.5)))
prior_prec_cell <- list(prec = list(prior = "loggamma", param = c(2, 0.0002)))

# fit linear model 
fit_inla_ask <- inla(resid ~ R + f(eqid, model = "iid", hyper = prior_prec_tau) + f(statid, model = "iid", hyper = prior_prec_phiS2S),
                     data = df.covar,
                     family="gaussian",
                     control.fixed = prior.fixed,
                     control.family = list(hyper = list(prec = prior_prec)),
                     control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config = TRUE))



# -----------------------------------------------
# spatial model for event and station terms terms
# -----------------------------------------------
# coordinates

X1_e <- coords_eq[eqid,2]
X2_e <- coords_eq[eqid,1]

X1_s <- coords_stat[statid,2]
X2_s <- coords_stat[statid,1]

coords <- unique(cbind(c(X1_e,X1_s), c(X2_e,X2_s)))

co_eq <- cbind(X1_e,X2_e)
co_eq2 <- unique(co_eq)

co_stat <- cbind(X1_s,X2_s)
co_stat2 <- unique(co_stat)

# build mesh
max.edge = 0.05
bound.outer = 0.5
mesh1 = inla.mesh.2d(loc=as.matrix(coords),
                     max.edge = c(1,5)*max.edge,
                     # - use 5 times max.edge in the outer extension/offset/boundary
                     cutoff = max.edge,
                     offset = c(5 * max.edge, bound.outer))

plot(mesh1)
points(co_stat[,1], co_stat[,2], col="red")
points(co_eq[,1], co_eq[,2], col="blue")
axis(1); axis(2)


# -----------------------------------------------------------------------
# define SPDE model

prior_psi_eq <- c(0.23,0.01)
prior_range_eq <- c(0.45,0.9)
# se prior and define model for spatally correlated event terms
spde_eq <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter
  mesh = mesh1, alpha = 2,
  # P(practic.range < 0.3) = 0.5
  prior.range = prior_range_eq,
  # P(sigma > 1) = 0.01
  prior.sigma = prior_psi_eq) 


prior_psi_stat <- c(0.23,0.01)
prior_range_stat <- c(0.45,0.9)
# se prior and define model for spatally correlated event terms
spde_stat <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter
  mesh = mesh1, alpha = 2,
  # P(practic.range < 0.3) = 0.5
  prior.range = prior_range_stat,
  # P(sigma > 1) = 0.01
  prior.sigma = prior_psi_stat)

prior_psi_vs <- c(0.23,0.01)
prior_range_vs <- c(0.45,0.9)
# se prior and define model for spatally correlated event terms
spde_vs <- inla.spde2.pcmatern(
  # Mesh and smoothness parameter
  mesh = mesh1, alpha = 2,
  # P(practic.range < 0.3) = 0.5
  prior.range = prior_range_vs,
  # P(sigma > 1) = 0.01
  prior.sigma = prior_psi_vs)

# make A matrices
A_eq <- inla.spde.make.A(mesh1, loc = co_eq)
idx.eq <- inla.spde.make.index("idx.eq",spde_eq$n.spde)

A_stat <- inla.spde.make.A(mesh1, loc = co_stat)
idx.stat <- inla.spde.make.index("idx.stat",spde_stat$n.spde)

A_vs <- inla.spde.make.A(mesh1, loc = co_stat)
idx.vs <- inla.spde.make.index("idx.vs",spde_vs$n.spde, weights = lnVS)


# create the stack
stk_eq_statvs_cell <- inla.stack(
  data = list(y = resid),
  A = list(A_eq, A_stat, A_vs, 1), 
  effects = list(idx.eq = idx.eq,
                 idx.stat = idx.stat,
                 idx.vs = idx.vs,
                 df.covar
  ),
  tag = 'model_eq_statvs_cell')

# define formula
form <- y ~ 0 + intercept + R +
  f(eq, model = "iid", hyper = prior_prec_tau) + f(stat, model = "iid", hyper = prior_prec_phiS2S) +
  f(idx.eq, model = spde_eq) + f(idx.stat, model = spde_stat) +
  f(idx.vs, model = spde_vs) +
  f(idx_cell, model = "z", Z = as.matrix(RC), hyper = prior_prec_cell)

# run spatial model
# use config=TRUE to be able to sample from the model later on
fit_inla_spatialvs_ask <- inla(form, 
                             data = inla.stack.data(stk_eq_statvs_cell), 
                             control.predictor = list(A = inla.stack.A(stk_eq_statvs_cell)),
                             control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config = TRUE),
                             family="gaussian",
                             control.fixed = prior.fixed,
                             control.family = list(hyper = list(prec = prior_prec)),
)
