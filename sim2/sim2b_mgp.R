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
a <- 10
c <- 20
kappa <- 1
sig_sq <- 150
nu <- 1.5
tau_sq <- 0.1

# assign partitions
n_easting <- 2
n_northing <- 6

# mcmc
mcmc <- list(save = 1000, burn = 5000, thin = 2)

# prior
tmaxdist <- max(dist(coords[,3], method = "manhattan"))
la <- (1/0.9-1)/tmaxdist                                                        # correlation = 0.9 at the maximum distance
ua <- (1/0.4-1)/(0.1*tmaxdist)                                                  # correlation = 0.4 at 1/10 of the maximal distance
spmaxdist <- max(dist(coords[,1:2]))
lc <- -log(0.95)/(0.1*spmaxdist)                                                # correlation = 0.95 at 1/10 of the maximal distance
uc <- -log(0.05)/(0.1*spmaxdist)                                                # correlation = 0.05 at 1/10 of the maximal distance

# generate seed
set.seed(123)
seedsave <- sample(10000, 100)
seedsave <- c(seedsave[51:74], 3)

# save results
rmspe = mape = coverage = meanwidth = tausq_hat = sigsq_hat <- rep(0, 25)
psi_hat = tot_time <- matrix(0, nrow = 3, ncol = 25)

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
  # mgp
  y <- data$y
  y[-tr_idx] <- NA
  x <- matrix(1, nrow = n, ncol = 1)

  time_mgp <- system.time({
    meshout <- spmeshed(family = c("gaussian"),
                        y = y, x = x,
                        coords = data[, c("easting", "northing", "time")],
                        k = 1,
                        axis_partition = c(n_easting, n_northing, n_time),
                        prior = list(phi = c(min(lc,la), max(uc,ua))),
                        settings = list(adapting = TRUE, forced_grid = NULL,
                                        cache = NULL, ps = TRUE, saving = TRUE),
                        starting = list(beta = NULL, tausq = NULL, theta = NULL,
                                        lambda = matrix(1,1,1), w = NULL,
                                        nu = NULL, mcmcsd = .05,
                                        mcmc_startfrom = 0),
                        debug = list(sample_beta = TRUE, sample_tausq = TRUE,
                                     sample_theta = TRUE, sample_w = TRUE,
                                     sample_lambda = FALSE,
                                     verbose = FALSE, debug = FALSE),
                        n_samples = mcmc$save,
                        n_burn = mcmc$burn,
                        n_thin = mcmc$thin,
                        verbose = 0,
                        n_threads = 10)
  })

  # total time spent in minutes
  tot_time[,s] <- (time_mgp/60)[1:3]

  # parameter estimation summary
  tausq_mcmc <- meshout$tausq_mcmc[,mcmc$thin*(1:mcmc$save)]
  psi_mcmc <- meshout$theta_mcmc[-4,,mcmc$thin*(1:mcmc$save)]                   # temporal decay, spatial decay, interaction
  sigsq_mcmc <- meshout$theta_mcmc[4,,mcmc$thin*(1:mcmc$save)]

  tausq_hat[s] <- mean(tausq_mcmc)
  psi_hat[,s] <- rowMeans(psi_mcmc)
  sigsq_hat[s] <- mean(sigsq_mcmc)

  # prediction summary
  y_mgp <- meshout$yhat_mcmc %>% summary_list_mean()
  y_mgp_low <- meshout$yhat_mcmc %>% summary_list_q(q = 0.025)
  y_mgp_hi <- meshout$yhat_mcmc %>% summary_list_q(q = 0.975)
  mgp_df_out <- data.frame(rbind(coords_tr[,c("easting", "northing", "time")],
                                 coords_tt[,c("easting", "northing", "time")]),
                           y = c(y_tr, y_tt)) %>%
    left_join(cbind(meshout$coordsdata, y_mgp, y_mgp_low, y_mgp_hi))

  mgp_res <- mgp_df_out[-c(1:n_tr),] %>%
    mutate(cover = ifelse((y > y_mgp_low & y < y_mgp_hi), 1, 0),
           range = y_mgp_hi - y_mgp_low) %>%
    summarise(rmspe = sqrt(mean((y-y_mgp)^2)),
              mape = mean(abs(y-y_mgp)),
              coverage = mean(cover),
              meanwidth = mean(range))

  coverage[s] <- mgp_res$coverage
  meanwidth[s] <- mgp_res$meanwidth
  rmspe[s] <- mgp_res$rmspe
  mape[s] <- mgp_res$mape

  setTxtProgressBar(pb, s/25)
}
close(pb)

saveRDS(list(rmspe = rmspe,
             mape = mape,
             coverage = coverage,
             meanwidth = meanwidth,
             tot_time = tot_time,
             tausq_hat = tausq_hat,
             psi_hat = psi_hat,
             sigsq_hat = sigsq_hat,
             meshout = meshout),                                                # seed = 3
        paste0(path, "sim2/sim2b_mgp_all.RDS"))

saveRDS(list(rmspe = rmspe,
             mape = mape,
             coverage = coverage,
             meanwidth = meanwidth),
        paste0(path, "sim2/sim2b_mgp.RDS"))
