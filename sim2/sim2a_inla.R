###########################
# Simulation: Northerlies #
###########################
#########
# Model #
#########
# y(l) = w(l) + eps(l), i = 1, ..., n
# l = (s, t)
# eps(l) ~ N(0, tau_sq)
# Assume two spatial parents (not fixed across times)
# Assume a nonseparable, stationary, and isotropic spatio-temporal Matern correlation function
# a: temporal decay, a > 0
# c: spatial decay, c > 0
# kappa: interaction between space and time
# nu: smoothness

rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(INLA)
library(meshed) # meshed_0.2.tar

# path
path <- "~/GBAGs/"

########
# Data #
########
# specify number of grid on each axis to generate data
ngrid0 <- 193
n_time0 <- 59
xgrid0 <- seq(0, 1, length = ngrid0)
tgrid0 <- seq(0, 1, length = n_time0)
coords0 <- expand.grid(easting = xgrid0, northing = xgrid0, time = tgrid0) %>%
  arrange(time, easting, northing)
n0 <- nrow(coords0)

# for a subset of data to fit
ngrid <- 25
n_time <- 30
xgrid <- seq(0, 1, length = ngrid)
tgrid <- seq(0, 1, length = n_time)
coords <- expand.grid(easting = xgrid, northing = 1-xgrid, time = tgrid) %>%
  arrange(time, easting, northing)
n <- nrow(coords)

# true parameter values
a <- 5
c <- 20
kappa <- 1
sig_sq <- 150
nu <- 1.5
tau_sq <- 0.1

# assign partitions
n_easting <- 2
n_northing <- 6

# generate seed
set.seed(123)
seedsave <- sample(10000, 100)
seedsave <- c(seedsave[51:74], 3)

# save results
## stationary
rmspe = mape = coverage = meanwidth = tausq_hat = sigsq_hat <- rep(0, 25)
tot_time <- matrix(0, 3, 25)

## nonstationary
rmspe_ns = mape_ns = coverage_ns = meanwidth_ns =
  tausq_hat_ns = sigsq_hat_ns <- rep(0, 25)
tot_time_ns <- matrix(0, 3, 25)

pb = txtProgressBar(style=3,width=50)
for (s in 1:25) {
  seed <- seedsave[s]
  set.seed(seed)

  # create w and y
  data0 <- rmeshedgp(coords = coords0 %>%
                       rename(Var1 = easting, Var2 = northing, Var3 = time),
                     theta = c(a, c, kappa, sig_sq, nu),
                     axis_partition = c(ngrid0, 1, n_time0),
                     n_threads = 10) %>%
    mutate(easting = Var2, northing = 1-Var1, time = Var3) %>%
    select(easting, northing, time, w)
  data0$y <- data0$w + sqrt(tau_sq)*rnorm(n0)
  data <- inner_join(coords, data0, by = c("easting", "northing", "time"))

  #################
  # Training data #
  #################
  # train data
  prob <- 0.8
  grid0 <- expand.grid(easting = xgrid, northing = 1-xgrid)

  idx_tr_tmp <- sort(sample(nrow(grid0), nrow(grid0)*prob))
  xygrid_tr <- list()
  for (tt in 1:n_time) {
    xygrid_tr[[tt]] <- data.frame(grid0[idx_tr_tmp,], time = tgrid[tt])
  }
  xygrid_tr <- do.call(rbind, xygrid_tr) %>% arrange(easting, northing, time)
  data_tr <- left_join(xygrid_tr,
                       data %>% mutate(idx = 1:n),
                       by = c("easting", "northing", "time"))
  tr_idx <- data_tr$idx
  n_tr <- length(tr_idx)
  n_tt <- n - n_tr
  w_tr <- data_tr$w
  w_tt <- data$w[-tr_idx]
  y_tr <- data_tr$y
  y_tt <- data$y[-tr_idx]

  coords_tr <- data_tr[, c("easting", "northing", "time")]
  coords_tt <- data[-tr_idx, c("easting", "northing", "time")]

  ###########
  # Methods #
  ###########
  # INLA
  # https://becarioprecario.bitbucket.io/spde-gitbook/ch-spacetime.html#sec:stcont
  # mesh
  diff_east <- diff(range(coords_tr$easting))
  max_edge <- diff_east/15
  uniq_coords_tr <- unique(coords_tr[,1:2])
  prmesh1 <- inla.mesh.create.helper(points = uniq_coords_tr,
                                     max.edge = c(2,10)*max_edge,
                                     cutoff = max_edge/1.5,
                                     offset = c(max_edge, diff_east/7))

  # define the SPDE (Penalized Complexity-priors)
  # stationary
  spde <- inla.spde2.pcmatern(mesh = prmesh1,
                              prior.range = c(0.1, 0.01),                       # P(range < 0.1) = 0.01
                              prior.sigma = c(1, 0.01))                         # P(sigma > 1) = 0.01

  # nonstationary
  # https://becarioprecario.bitbucket.io/spde-gitbook/ch-nonstationarity.html
  # range vary by coordinates of mesh nodes
  me_max <- max(prmesh1$loc[,1])
  logkappa0 <- .5*log(8)
  logtau0 <- -.5*log(4*pi)
  logtau0 <- logtau0 - logkappa0
  spde_ns <- inla.spde2.matern(prmesh1,
                               B.tau = cbind(logtau0, -1, 1,
                                             prmesh1$loc[,1], prmesh1$loc[,2],
                                             (prmesh1$loc[,1]-1)*(prmesh1$loc[,2]-1)),
                               B.kappa = cbind(logkappa0, 0, -1,
                                               -prmesh1$loc[,1], -prmesh1$loc[,2],
                                               -(prmesh1$loc[,1]-1)*(prmesh1$loc[,2]-1)),
                               theta.prior.mean = rep(0, 5),
                               theta.prior.prec = rep(.3, 5))

  # for space-TIME model
  iset <- inla.spde.make.index('i', n.spde = spde$n.spde, n.group = n_time)
  iset_ns <- inla.spde.make.index('i', n.spde = spde_ns$n.spde, n.group = n_time)

  # connect observations to mesh nodes
  A.obs <- inla.spde.make.A(mesh = prmesh1,
                            loc = as.matrix(coords_tr[,1:2]),
                            group = as.numeric(
                              as.character(coords_tr$time*(n_time-1) + 1)))     # group should be an index
  ybar <- mean(y_tr)
  stk.obs <- inla.stack(
    data = list(y = y_tr - ybar),
    A = list(A.obs),
    effects = list(iset),
    tag = 'obs')
  stk.obs_ns <- inla.stack(
    data = list(y = y_tr - ybar),
    A = list(A.obs),
    effects = list(iset_ns),
    tag = 'obs')

  # same for prediction
  proj.pred <- inla.mesh.projector(prmesh1, loc = as.matrix(coords_tt[,1:2]))
  A.pred <- inla.spde.make.A(prmesh1, loc = proj.pred$loc,
                             group = as.numeric(
                               as.character(coords_tt$time*(n_time-1) + 1)))
  stk.pred <- inla.stack(data = list(y = NA),
                         A = list(A.pred),
                         effects = list(iset),
                         tag = "pred")
  stk.pred_ns <- inla.stack(data = list(y = NA),
                            A = list(A.pred),
                            effects = list(iset_ns),
                            tag = "pred")
  stk <- inla.stack(stk.obs, stk.pred)
  stk_ns <- inla.stack(stk.obs_ns, stk.pred_ns)

  # model formula
  h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0.05, 0.99)))          # P(temporal auto-corr > 0.05) = 0.99
  formulae <- y ~ -1 + f(i, model = spde, group = i.group,
                         control.group = list(model = 'ar1', hyper = h.spec))
  formulae_ns <- y ~ -1 + f(i, model = spde_ns, group = i.group,
                            control.group = list(model = 'ar1', hyper = h.spec))

  # PC prior on the autoregressive parameter
  prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))                     # P(precision of AR > 1) = 0.01

  # model fitting
  time_inla <- system.time({
    inlaout <- inla(formulae, data = inla.stack.data(stk),
                    control.predictor = list(compute = TRUE,
                                             A = inla.stack.A(stk)),
                    control.family = list(hyper = list(prec = prec.prior)),
                    num.threads = 10,
                    verbose = FALSE)
  })

  time_inla_ns <- system.time({
    inlaout_ns <- inla(formulae_ns, data = inla.stack.data(stk_ns),
                       control.predictor = list(compute = TRUE,
                                                A = inla.stack.A(stk_ns)),
                       control.family = list(hyper = list(prec = prec.prior)),
                       num.threads = 10,
                       verbose = FALSE)
  })

  # total time spent in minutes
  tot_time[,s] <- (time_inla/60)[1:3]
  tot_time_ns[,s] <- (time_inla_ns/60)[1:3]

  # parameter estimation summary
  ## stationary
  tausq_hat[s] <- mean(1/inlaout$marginals.hyperpar$
                         `Precision for the Gaussian observations`[,1])

  ## nonstationary
  tausq_hat_ns[s] <- mean(1/inlaout_ns$marginals.hyperpar$
                            `Precision for the Gaussian observations`[,1])

  # prediction summaries
  index.pred <- c(inla.stack.index(stk, "pred")$data)
  ## stationary
  y_inla <- ybar + inlaout$summary.fitted.values[index.pred, "mean"]
  inla_sd <- sqrt(tausq_hat[s])
  y_inla_low <- y_inla + qnorm(.025)*inla_sd
  y_inla_hi <- y_inla + qnorm(.975)*inla_sd

  coverage[s] <- mean(ifelse(y_tt > y_inla_low & y_tt < y_inla_hi, 1, 0))
  meanwidth[s] <- mean(y_inla_hi - y_inla_low)
  rmspe[s] <- sqrt(mean((y_tt-y_inla)^2))
  mape[s] <- mean(abs(y_tt-y_inla))

  ## nonstationary
  y_inla_ns <- ybar + inlaout_ns$summary.fitted.values[index.pred, "mean"]
  inla_sd_ns <- sqrt(tausq_hat_ns[s])
  y_inla_low_ns <- y_inla_ns + qnorm(.025)*inla_sd_ns
  y_inla_hi_ns <- y_inla_ns + qnorm(.975)*inla_sd_ns

  coverage_ns[s] <- mean(ifelse(y_tt > y_inla_low_ns & y_tt < y_inla_hi_ns, 1, 0))
  meanwidth_ns[s] <- mean(y_inla_hi_ns - y_inla_low_ns)
  rmspe_ns[s] <- sqrt(mean((y_tt-y_inla_ns)^2))
  mape_ns[s] <- mean(abs(y_tt-y_inla_ns))

  setTxtProgressBar(pb, s/25)
}
close(pb)

saveRDS(list(rmspe = rmspe,
             mape = mape,
             coverage = coverage,
             meanwidth = meanwidth,
             tot_time = tot_time,
             tausq_hat = tausq_hat,
             inlaout = inlaout,
             rmspe_ns = rmspe_ns,
             mape_ns = mape_ns,
             coverage_ns = coverage_ns,
             meanwidth_ns = meanwidth_ns,
             tot_time_ns = tot_time_ns,
             tausq_hat_ns = tausq_hat_ns,
             inlaout_ns = inlaout_ns,
             data = list(prmesh = prmesh1,
                         iset = iset,
                         iset_ns = iset_ns)),                                   # seed = 3
        paste0(path, "sim2/sim2a_inla_all.RDS"))

saveRDS(list(rmspe = rmspe,
             mape = mape,
             coverage = coverage,
             meanwidth = meanwidth,
             rmspe_ns = rmspe_ns,
             mape_ns = mape_ns,
             coverage_ns = coverage_ns,
             meanwidth_ns = meanwidth_ns),
        paste0(path, "sim2/sim2a_inla.RDS"))
