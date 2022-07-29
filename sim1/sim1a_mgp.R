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

# mcmc
mcmc <- list(save = 1000, burn = 10000, thin = 7)

# true parameters
tau_sq <- 0.01
sig_sq <- 2
a <- 5
c <- 0.5
kappa <- 0.9
beta <- 2

# prior
tmaxdist <- max(dist(coords[,3], method = "manhattan"))
la <- (1/0.2-1)/tmaxdist                                                        # correlation = 0.2 at the maximum distance
ua <- (1/0.2-1)/(0.5*tmaxdist)                                                  # correlation = 0.2 at the half maximal distance
spmaxdist <- max(dist(coords[,1:2]))
lc <- -log(0.8)/(spmaxdist)                                                     # correlation = 0.8 at the maximum distance
uc <- -log(0.8)/(0.2*spmaxdist)                                                 # correlation = 0.8 at one fifth of the maximum distance

# generate seed
set.seed(123)
seedsave <- sample(10000, 100)
seedsave <- c(seedsave[1:24], 123)

# save results
rmspe = mape = coverage = meanwidth <- rep(0, 25)
beta_hat = beta_ci = beta_ess =
  tausq_hat = tausq_ci = tausq_ess =
  sigsq_hat = sigsq_ci = sigsq_ess = sigsqc_ess <- rep(0, 25)
psi_hat = psi_ci = psi_ess <- matrix(0, nrow = 3, ncol = 25)
time_per_iter <- rep(0, 25)
y_pred_ess <- list()

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
  # mgp - partitions do not match
  y <- data$y
  y[-tr_idx] <- NA
  X <- as.matrix(data$X)

  meshout <- spmeshed(family = c("gaussian"),
                      y = y, x = X, coords = coords, k = 1,
                      axis_partition = c(n_easting, n_northing, n_time),
                      grid_custom = list(axis_interval_partition =
                                           list(easting = easting_midval,
                                                northing = northing_midval,
                                                time = time_midval)),
                      prior = list(phi = c(min(lc,la), max(uc,ua))),
                      settings = list(adapting = TRUE, forced_grid = NULL,
                                      cache = NULL, ps = TRUE, saving = TRUE),
                      starting = list(beta = NULL, tausq = NULL, theta = NULL,
                                      lambda = matrix(1,1,1), w = NULL,
                                      nu = NULL, mcmcsd = .05, mcmc_startfrom = 0),
                      debug = list(sample_beta = TRUE, sample_tausq = TRUE,
                                   sample_theta = TRUE, sample_w = TRUE,
                                   sample_lambda = FALSE,
                                   verbose = FALSE, debug = FALSE),
                      n_samples = mcmc$save,
                      n_burn = mcmc$burn,
                      n_thin = mcmc$thin,
                      verbose = 0,
                      n_threads = 10)

  # time spent per iteration in seconds
  time_per_iter[s] <- meshout$mcmc_time/meshout$mcmc

  # parameter estimation summary
  beta_mcmc <- meshout$beta_mcmc[,,mcmc$thin*(1:mcmc$save)]
  tausq_mcmc <- meshout$tausq_mcmc[,mcmc$thin*(1:mcmc$save)]
  psi_mcmc <- meshout$theta_mcmc[-4,,mcmc$thin*(1:mcmc$save)]                   # temporal decay, spatial decay, interaction
  sigsq_mcmc <- meshout$theta_mcmc[4,,mcmc$thin*(1:mcmc$save)]

  beta_hat[s] <- mean(beta_mcmc)
  tausq_hat[s] <- mean(tausq_mcmc)
  psi_hat[,s] <- rowMeans(psi_mcmc)
  sigsq_hat[s] <- mean(sigsq_mcmc)

  # 95% CI coverage
  beta_ci[s] <- (beta > quantile(beta_mcmc, 0.025) &
                   beta < quantile(beta_mcmc, 0.975))
  tausq_ci[s] <- (tau_sq > quantile(tausq_mcmc, 0.025) &
                    tau_sq < quantile(tausq_mcmc, 0.975))
  psi_ci[1,s] <- (a > quantile(psi_mcmc[1,], 0.025) &
                    a < quantile(psi_mcmc[1,], 0.975))
  psi_ci[2,s] <- (c > quantile(psi_mcmc[2,], 0.025) &
                    c < quantile(psi_mcmc[2,], 0.975))
  psi_ci[3,s] <- (kappa > quantile(psi_mcmc[3,], 0.025) &
                    kappa < quantile(psi_mcmc[3,], 0.975))
  sigsq_ci[s] <- (sig_sq > quantile(sigsq_mcmc, 0.025) &
                    sig_sq < quantile(sigsq_mcmc, 0.975))

  # effective sample size
  beta_ess[s] <- coda::effectiveSize(beta_mcmc)
  tausq_ess[s] <- coda::effectiveSize(tausq_mcmc)
  sigsq_ess[s] <- coda::effectiveSize(sigsq_mcmc)
  psi_ess[,s] <- apply(psi_mcmc, 1, coda::effectiveSize)
  sigsqc_ess[s] <- coda::effectiveSize(sigsq_mcmc*psi_mcmc[2,])
  mgp_df_out <- data.frame(rbind(coords_tr[,c("easting", "northing", "time")],
                                 coords_tt[,c("easting", "northing", "time")])) %>%
    left_join(cbind(meshout$coordsdata, do.call(cbind, meshout$yhat_mcmc)),
              by = c("easting", "northing", "time"))
  y_pred_ess[[s]] <- apply(mgp_df_out[-c(1:n_tr), -c(1:4)],
                           1, coda::effectiveSize) %>% as.numeric()

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
             time_per_iter = time_per_iter,
             beta_hat = beta_hat,
             tausq_hat = tausq_hat,
             psi_hat = psi_hat,
             sigsq_hat = sigsq_hat,
             beta_ci = beta_ci,
             tausq_ci = tausq_ci,
             psi_ci = psi_ci,
             sigsq_ci = sigsq_ci,
             beta_ess = beta_ess,
             tausq_ess = tausq_ess,
             psi_ess = psi_ess,
             sigsq_ess = sigsq_ess,
             sigsqc_ess = sigsqc_ess,
             y_pred_ess = y_pred_ess,
             meshout = meshout),                                                # seed = 123
        paste0(path, "sim1/sim1a_mgp_all.RDS"))

saveRDS(list(rmspe = rmspe,
             mape = mape,
             coverage = coverage,
             meanwidth = meanwidth,
             time_per_iter = time_per_iter,
             beta_hat = beta_hat,
             tausq_hat = tausq_hat,
             psi_hat = psi_hat,
             sigsq_hat = sigsq_hat,
             beta_ci = beta_ci,
             tausq_ci = tausq_ci,
             psi_ci = psi_ci,
             sigsq_ci = sigsq_ci,
             beta_ess = beta_ess,
             tausq_ess = tausq_ess,
             psi_ess = psi_ess,
             sigsq_ess = sigsq_ess,
             sigsqc_ess = sigsqc_ess,
             y_pred_ess = y_pred_ess),
        paste0(path, "sim1/sim1a_mgp.RDS"))


