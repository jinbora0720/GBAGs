# Bag of DAGs Simulations #
#########
# Model #
#########
# y(l) = X(l)^T %*% beta + w(l) + eps(l), i = 1, ..., n
# l = (s, t)
# X(l)^T: (1 x p),
# eps(l) ~ N(0, tau_sq)
# w|z=h ~ N(0, sig_sq*Ctilde_h(phi)) for h = 1, ..., K where
# Ctilde_h(phi)^{-1} = (I-H_h(phi))^T %*% R_h(phi)^{-1} %*% (I-H_h(phi))
# H_h(phi) is a sparse matrix whose non-zero elements' locations are determined by the hth DAG.
# R_h(phi) is a diagonal matrix.
# H_h(phi) and R_h(phi) are conditional mean coefficient and conditional variance of normal distributions.
# Assume a nonseparable, stationary, and isotropic spatio-temporal correlation function
# Cor(w(l), w(l + (h,u))) = exp(-c||h||/(a|u|+1)^{kappa/2})/(a|u|+1)
# c: spatial decay, c > 0
# a: temporal decay, a > 0
# kappa: interaction between space and time

rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(bags)
library(meshed)
library(INLA)

# path
path <- "~/GBAGs/"

########
# Data #
########
# specify number of grid on each axis
ngrid <- 40
n_time <- 8
xgrid <- seq(0, 1, length = ngrid)
ygrid <- xgrid
tgrid <- seq(0, 1, length = n_time)
xytgrid <- expand.grid(easting = xgrid, northing = ygrid, time = tgrid) %>%
  arrange(time, easting, northing)
n <- nrow(xytgrid)
coords <- xytgrid

# assign partitions (irregular partitions)
n_easting <- 6
n_northing <- 6

# directions
directions <- c("W", "NW", "N", "NE")

# data
nd <- floor(log10(max(n_northing, n_easting, n_time))) + 1
format <- paste0("%0", nd, ".0f,","%0", nd, ".0f,", "%0", nd, ".0f")

coords_ptt <- coords
northing_sint <- seq(0, 1, length = (n_northing-2)^2+1)
northing_Lint <- seq(0, 1, length = (n_northing-2)+1)
breaks_northing <- c(northing_sint[1:2], northing_Lint[2:(n_northing-2)],
                     northing_sint[(n_northing-2)^2+0:1])
coords_ptt$row <- (n_northing+1) -
  as.numeric(cut(coords_ptt$northing,
                 breaks = breaks_northing,
                 labels = 1:n_northing, include.lowest = T)) # 1 from the top
easting_sint <- seq(0, 1, length = (n_easting-2)^2+1)
easting_Lint <- seq(0, 1, length = (n_easting-2)+1)
breaks_easting <- c(easting_sint[1:2], easting_Lint[2:(n_easting-2)],
                    easting_sint[(n_easting-2)^2+0:1])
coords_ptt$col <- as.numeric(cut(coords_ptt$easting,
                                 breaks = breaks_easting,
                                 labels = 1:n_easting, include.lowest = T))
coords_ptt$time_d <- as.numeric(as.character(coords_ptt$time*(n_time-1) + 1))
coords_ptt$partition <- sprintf(format,
                                coords_ptt$row, coords_ptt$col,
                                coords_ptt$time_d)

ptts <- sort(unique(coords_ptt$partition))
z_true <- rep("W", length(ptts))
ptt_NW <- grep("\\d+,4,2|\\d+,5,2|\\d+,6,2|\\d+,\\d+,3|1,\\d+,4|2,\\d+,4|3,\\d+,4|4,2,8|4,3,8|4,4,8|4,5,8|4,6,8|5,6,8|6,6,8", ptts)
z_true[ptt_NW] <- "NW"
ptt_N <- grep("4,\\d+,4|5,\\d+,4|6,\\d+,4|\\d+,\\d+,5|1,\\d+,6|2,\\d+,6|3,\\d+,6|3,1,8|3,2,8|3,3,8|3,4,8|3,5,8|4,1,8", ptts)
z_true[ptt_N] <- "N"
ptt_NE <- grep("4,\\d+,6|5,\\d+,6|6,\\d+,6|\\d+,\\d+,7|5,1,8|5,2,8|5,3,8|5,4,8|5,5,8|6,1,8|6,2,8|6,3,8|6,4,8|6,5,8", ptts)
z_true[ptt_NE] <- "NE"

# true parameters
tau_sq <- 0.01
sig_sq <- 2
a <- 10
c <- 0.1
kappa <- 0.2
beta <- 2

# generate seed
set.seed(123)
seedsave <- sample(10000, 100)
seedsave <- seedsave[26:50]

# save results
## stationary cov
rmspe = mape = coverage = meanwidth <- rep(0, 25)
beta_hat = beta_ci =
  tausq_hat = tausq_ci <- rep(0, 25)
tot_time <- matrix(0, nrow = 3, ncol = 25)

## nonstationary
rmspe_ns = mape_ns = coverage_ns = meanwidth_ns <- rep(0, 25)
beta_hat_ns = beta_ci_ns =
  tausq_hat_ns = tausq_ci_ns <- rep(0, 25)
tot_time_ns <- matrix(0, nrow = 3, ncol = 25)

pb = txtProgressBar(style=3,width=50)
for (s in 1:25) {
  seed <- seedsave[s]
  set.seed(seed)

  data <- rbag(coords = as.matrix(coords),
               n_partition = c(n_easting, n_northing, n_time),
               breaks_partition = list(breaks_easting = breaks_easting,
                                       breaks_northing = breaks_northing,
                                       breaks_time = NULL),
               z = z_true,
               params = list(tau_sq = tau_sq,
                             sig_sq = sig_sq,
                             a = a, c = c, kappa = kappa,
                             beta = beta),
               seed = seed)

  #################
  # Training data #
  #################
  # train data
  prob <- 0.8
  grid0 <- expand.grid(easting = xgrid, northing = ygrid)

  idx_tr_tmp <- sort(sample(nrow(grid0), nrow(grid0)*prob))
  xygrid_tr <- list()
  for (tt in 1:n_time) {
    xygrid_tr[[tt]] <- data.frame(grid0[idx_tr_tmp,], time = tgrid[[tt]])
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
  X_tr <- data_tr$X
  X_tt <- data$X[-tr_idx]

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
                              prior.range = c(9, 0.01),                         # P(range = 10 < 9) = 0.01
                              prior.sigma = c(1.5, 0.01))                       # P(sigma = 1.2 > 1.5) = 0.01

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
                            group = as.numeric(as.character(
                              coords_tr$time*(n_time-1) + 1)))                  # group should be an index
  ybar <- mean(y_tr)
  stk.obs <- inla.stack(
    data = list(y = y_tr - ybar),
    A = list(A.obs, 1),
    effects = list(iset, X = X_tr),
    tag = 'obs')
  stk.obs_ns <- inla.stack(
    data = list(y = y_tr - ybar),
    A = list(A.obs, 1),
    effects = list(iset_ns, X = X_tr),
    tag = 'obs')

  # same for prediction
  proj.pred <- inla.mesh.projector(prmesh1, loc = as.matrix(coords_tt[,1:2]))
  A.pred <- inla.spde.make.A(prmesh1, loc = proj.pred$loc,
                             group = as.numeric(as.character(
                               coords_tt$time*(n_time-1) + 1)))
  stk.pred <- inla.stack(data = list(y = NA),
                         A = list(A.pred, 1),
                         effects = list(iset, X = X_tt),
                         tag = "pred")
  stk.pred_ns <- inla.stack(data = list(y = NA),
                            A = list(A.pred, 1),
                            effects = list(iset_ns, X = X_tt),
                            tag = "pred")
  stk <- inla.stack(stk.obs, stk.pred)
  stk_ns <- inla.stack(stk.obs_ns, stk.pred_ns)

  # model formula
  h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0.05, 0.99)))          # P(temporal auto-corr > 0.05) = 0.99
  formulae <- y ~ -1 + X + f(i, model = spde, group = i.group,
                             control.group = list(model = 'ar1', hyper = h.spec))
  formulae_ns <- y ~ -1 + X + f(i, model = spde_ns, group = i.group,
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
  beta_hat[s] <- as.numeric(inlaout$summary.fixed["mean"])
  tausq_summary <- inlaout$internal.marginals.hyperpar[1] %>%                   # log precision for the Gaussian observations
    lapply(function(m) {
      inla.tmarginal(function(x) 1/exp(x), m)
    }) %>%
    sapply(function(m)
      unlist(inla.zmarginal(m, silent = TRUE)))
  tausq_hat[s] <- tausq_summary["mean",]

  ## nonstationary
  beta_hat_ns[s] <- as.numeric(inlaout_ns$summary.fixed["mean"])
  tausq_summary_ns <- inlaout_ns$internal.marginals.hyperpar[1] %>%             # log precision for the Gaussian observations
    lapply(function(m) {
      inla.tmarginal(function(x) 1/exp(x), m)
    }) %>%
    sapply(function(m)
      unlist(inla.zmarginal(m, silent = TRUE)))
  tausq_hat_ns[s] <- tausq_summary_ns["mean",]

  # 95% CI coverage
  ## stationary
  beta_ci[s] <- (beta > inlaout$summary.fixed["0.025quant"] &
                   beta < inlaout$summary.fixed["0.975quant"]) %>% as.numeric()
  tausq_ci[s] <- (tau_sq > tausq_summary["quant0.025",] &
                    tau_sq < tausq_summary["quant0.975",])

  ## nonstationary
  beta_ci_ns[s] <- (beta > inlaout_ns$summary.fixed["0.025quant"] &
                      beta < inlaout_ns$summary.fixed["0.975quant"]) %>%
    as.numeric()
  tausq_ci_ns[s] <- (tau_sq > tausq_summary_ns["quant0.025",] &
                       tau_sq < tausq_summary_ns["quant0.975",])

  # prediction summaries
  index.pred <- c(inla.stack.index(stk, "pred")$data)
  ## stationary
  y_inla <- ybar + inlaout$summary.fitted.values[index.pred, "mean"]            # include Xb
  inla_sd <- sqrt(inlaout$summary.fitted.values[index.pred, "sd"]^2 +           # variability from Xb fitting
                    tausq_hat[s])
  y_inla_low <- y_inla + qnorm(.025)*inla_sd
  y_inla_hi <- y_inla + qnorm(.975)*inla_sd

  coverage[s] <- mean(ifelse(y_tt > y_inla_low & y_tt < y_inla_hi, 1, 0))
  meanwidth[s] <- mean(y_inla_hi - y_inla_low)
  rmspe[s] <- sqrt(mean((y_tt-y_inla)^2))
  mape[s] <- mean(abs(y_tt-y_inla))

  ## nonstationary
  y_inla_ns <- ybar + inlaout_ns$summary.fitted.values[index.pred, "mean"]
  inla_sd_ns <- sqrt(inlaout_ns$summary.fitted.values[index.pred, "sd"]^2 +
                       tausq_hat_ns[s])
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
             beta_hat = beta_hat,
             beta_ci = beta_ci,
             tausq_hat = tausq_hat,
             tausq_ci = tausq_ci,
             inlaout = inlaout,
             rmspe_ns = rmspe_ns,
             mape_ns = mape_ns,
             coverage_ns = coverage_ns,
             meanwidth_ns = meanwidth_ns,
             tot_time_ns = tot_time_ns,
             beta_hat_ns = beta_hat_ns,
             beta_ci_ns = beta_ci_ns,
             tausq_hat_ns = tausq_hat_ns,
             tausq_ci_ns = tausq_ci_ns,
             inlaout_ns = inlaout_ns),
        paste0(path, "sim1/sim1b_inla.RDS"))

