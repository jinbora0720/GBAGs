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
library(bags)

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

# mcmc
mcmc <- list(save = 1000, burn = 5000, thin = 2)

# z Chisq test
prop <- 0.35
n_samp <- ceiling(0.25*mcmc$save)
n1 = n2 <- n_samp

# prior
tmaxdist <- max(dist(coords[,3], method = "manhattan"))
la <- (1/0.9-1)/tmaxdist                                                        # correlation = 0.9 at the maximum distance
ua <- (1/0.1-1)/tmaxdist                                                        # correlation = 0.1 at the maximum distance
spmaxdist <- max(dist(coords[,1:2]))
lc <- -log(0.95)/(0.1*spmaxdist)                                                # correlation = 0.95 at 1/10 of the maximumal distance
uc <- -log(0.05)/(0.1*spmaxdist)                                                # correlation = 0.05 at 1/10 of the maximumal distance

# directions
directions <- c("NW", "N", "NE", "E")

# generate seed
set.seed(123)
seedsave <- sample(10000, 100)
seedsave <- c(seedsave[51:74], 3)

# save results
rmspe = mape = coverage = meanwidth = tausq_hat = sigsq_hat <- rep(0, 25)
psi_hat <- matrix(0, nrow = 3, ncol = 25)
est_time_per_iter = pred_time_per_iter <- matrix(0, nrow = 3, ncol = 25)
rownames(est_time_per_iter) = rownames(pred_time_per_iter) <-
  c("user", "system", "elapsed")
z_pvalue = z_save = y_pred_ess <- list()

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
  # bag
  save_data <- FALSE
  if (s == 25) {
    save_data <- TRUE
  }
  out <- bag(y = y_tr, X = NULL,
             coords = as.matrix(coords_tr[,c("easting", "northing", "time")]),
             X_pred = NULL,
             coords_pred = as.matrix(coords_tt[,c("easting", "northing", "time")]),
             n_partition = c(n_easting, n_northing, n_time),
             breaks_partition = list(breaks_easting = NULL,
                                     breaks_northing = NULL,
                                     breaks_time = NULL),
             directions = directions,
             init = list(tau_sq = NULL,
                         sig_sq = NULL,
                         w = NULL,
                         z = NULL,
                         psi = NULL,
                         Sn = NULL),
             hyper = list(at = NULL, bt = NULL,
                          as = NULL, bs = NULL,
                          la = la, ua = ua,
                          lc = lc, uc = uc,
                          mu0 = NULL, invV0 = NULL),
             mcmc = mcmc,
             n_threads = 10,
             seed = seed,
             verbose = FALSE,
             save_data = save_data,
             save_est = FALSE,
             debug = list(psi_fixed = FALSE, z_fixed = FALSE))

  # time spent per iteration
  est_time_per_iter[,s] <- as.numeric(out$est_time/out$est_iter)[1:3]
  pred_time_per_iter[,s] <- as.numeric(out$pred_time/out$pred_iter)[1:3]

  # parameter estimation summary
  tausq_hat[s] <- mean(out$tau_sq_save)
  psi_hat[,s] <- rowMeans(out$psi_save)
  sigsq_hat[s] <- mean(out$sig_sq_save)
  z_save[[s]] <- out$z_save

  # effective sample size
  y_pred_ess[[s]] <- apply(out$y_pred_save, 1, coda::effectiveSize)
  pvalue <- rep(0, nrow(out$z_save))
  for (ii in 1:nrow(out$z_save)) {
    samp1 <- out$z_save[ii,1:n_samp]
    samp2 <- out$z_save[ii,(mcmc$save - n_samp) + 1:n_samp]

    N1 <- c(sum(samp1 == directions[1]), sum(samp1 == directions[2]),
            sum(samp1 == directions[3]), sum(samp1 == directions[4]))
    p1 <- N1/n1
    N2 <- c(sum(samp2 == directions[1]), sum(samp2 == directions[2]),
            sum(samp2 == directions[3]), sum(samp2 == directions[4]))
    p2 <- N2/n2
    p <- (N1 + N2)/(n1 + n2)
    R <- which(p != 0)
    nR <- sum(p != 0)
    X_sq <- n1*sum((p1[R] - p[R])^2/p[R]) + n2*sum((p2[R] - p[R])^2/p[R])
    pvalue[ii] <- pchisq(X_sq, df = (nR-1)*(2-1), lower.tail = FALSE)
  }
  z_pvalue[[s]] <- pvalue

  # prediction summary
  y_bag <- rowMeans(out$y_pred_save)
  y_bag_qt <- apply(out$y_pred_save, 1, function(x)
    quantile(x, probs = c(0.025, 0.975)))
  y_bag_low <- y_bag_qt[1,]
  y_bag_hi <- y_bag_qt[2,]

  coverage[s] <- mean(ifelse(y_tt > y_bag_low & y_tt < y_bag_hi, 1, 0))
  meanwidth[s] <- mean(y_bag_hi - y_bag_low)
  rmspe[s] <- sqrt(mean((y_tt-y_bag)^2))
  mape[s] <- mean(abs(y_tt-y_bag))

  setTxtProgressBar(pb, s/25)
}
close(pb)

saveRDS(list(rmspe = rmspe,
             mape = mape,
             coverage = coverage,
             meanwidth = meanwidth,
             est_time_per_iter = est_time_per_iter,
             pred_time_per_iter = pred_time_per_iter,
             tausq_hat = tausq_hat,
             psi_hat = psi_hat,
             sigsq_hat = sigsq_hat,
             z_save = z_save,
             z_pvalue = z_pvalue,
             y_pred_ess = y_pred_ess,
             data0 = data0,                                                     # seed = 3
             out = out),
        paste0(path, "sim2/sim2a_bag.RDS"))

